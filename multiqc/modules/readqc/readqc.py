#!/usr/bin/env python

"""MultiQC module to parse output from ReadQC"""

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config
from multiqc.plots import table

import logging
import collections
import xml.etree.cElementTree
import re

# Initialise the logger
log = logging.getLogger(__name__)


class QcmlMultiqcModule(BaseMultiqcModule):

    def parse_qcml_by(self,qcml_contents,tag):
        """Parse a qcML file and return key-value pairs from the quality parameter entries."""
        root = xml.etree.cElementTree.fromstring(qcml_contents)
        parameters = dict()

        for qp in root.findall(".//{http://www.prime-xs.eu/ms/qcml}%s"%tag):
            # skip n/a values
            if qp.attrib['value'].startswith('n/a'):
                continue

            # replace 'percentage' with '%'
            qp_name = re.sub(r' percentage$', ' %', qp.attrib['name'])

            try:
                parameters[qp_name] = float(qp.attrib['value'])
            except ValueError:
                parameters[qp_name] = qp.attrib['value']

            # add description and accession number of the parameter to the header
            self.qcml[qp_name] = dict()
            self.add_attribute(qp_name,qp,['description', 'accession'])
        return parameters

    def add_attribute(self, parameter_name,content,attributes):
        for attrib in attributes:
            if attrib in content.attrib:
                self.qcml[parameter_name][attrib] = content.attrib[attrib]
            else:
                self.qcml[parameter_name][attrib] = ''

    def get_read_type(self,qcml_contents):
        root = xml.etree.cElementTree.fromstring(qcml_contents)
        source_files = list()

        for qp in root.findall(".//{http://www.prime-xs.eu/ms/qcml}metaDataParameter"):
            # skip n/a values
            if qp.attrib['value'].startswith('n/a'):
                continue
            # collect source files
            if qp.attrib['name'] == 'source file':
                source_files.append(qp.attrib['value'])

        # check the existence of R1 and R2 reads in the source files
        r1_exist = False
        r2_exist = False
        for source_file in source_files:
            if "R1" in source_file:
                r1_exist = True
            elif "R2" in source_file:
                r2_exist = True
        if r1_exist and r2_exist:
            return "paired-end"
        elif r1_exist or r2_exist:
            return "single"
        else:
            return "unknown"

    def make_description(self, keynames):
        """"Create description string from qcML quality parameter key name."""
        if len(keynames) == 1:
            desc = "{:s}".format(self.qcml[keynames[0]]['description'])
        else:
            desc = "<ul>" + "".join(
                ["<li>{:s}</li>".format(self.qcml[key]['description']) for key in keynames]) + "</ul>"
        return desc

    @staticmethod
    def dict_ordered_subset(d, ks):
        """Return subset of a dictionary as an OrderedDict object, ignoring non-existent keys."""
        od = collections.OrderedDict()
        for k in ks:
            try:
                od[k] = d[k]
            except KeyError:
                pass
        return od


class MultiqcModule(QcmlMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='ReadQC',
                                            anchor='readqc',
                                            href="https://github.com/imgag/ngs-bits",
                                            info="calculates QC metrics on unprocessed NGS reads. Paired-end reads count as two reads, and single-end reads count as one read.")

        # quality parameters from qcML with name, accession, description
        self.qcml = dict()
        # qc data for each sample
        self.qcdata = dict()
        self.read_types = dict()
        # parse qcml files
        for f in self.find_log_files('readqc', filecontents=True, filehandles=False):
            self.add_data_source(f)
            s_name = self.clean_s_name(f['s_name'], f['root'])
            self.qcdata[s_name] = self.parse_qcml_by(f['f'], "qualityParameter")
            self.read_types[s_name] = self.get_read_type(f['f'])
        # ignore samples if requested
        self.qcdata = self.ignore_samples(self.qcdata)

        # warn if no samples found
        if len(self.qcdata) == 0:
            raise UserWarning

        # add bases sequenced key, derived from bases sequenced (MB)
        self.qcml.pop('bases sequenced (MB)')
        self.qcml['bases sequenced'] = dict()
        self.qcml['bases sequenced']['description'] = 'Bases sequenced in total.'

        # add cluster count key, derived from read count
        self.qcml['cluster count'] = dict()
        self.qcml['cluster count']['description'] = 'Total number of clusters.'

        for s_name, kv in self.qcdata.items():
            kv['bases sequenced'] = kv['bases sequenced (MB)'] * 1e6

            # calculate cluster count according to the read type
            if self.read_types[s_name] == "single":
                kv['cluster count'] = kv['read count']
            elif self.read_types[s_name] == "paired-end":
                kv['cluster count'] = kv['read count']/2

            kv.pop('bases sequenced (MB)')

        # prepare table headers, use name and description from qcML
        headers = {qp_key: {
            'namespace': "ReadQC",
            'title': qp_key,
            'description': qp_entry['description'],
        } for qp_key, qp_entry in self.qcml.items()}
        headers['read count'].update({'suffix': config.read_count_prefix, 'format': '{:,.2f}',
                                      'modify': lambda x: x * config.read_count_multiplier,
                                      'scale': 'Purples', 'placement': 10})
        headers['cluster count'].update({'suffix': config.cluster_count_prefix, 'format': '{:,.2f}',
                                         'modify': lambda x: x * config.cluster_count_multiplier,
                                         'scale': 'Purples', 'placement': 20})
        headers['bases sequenced'].update({'suffix': config.base_count_prefix, 'format': '{:,.2f}',
                                           'modify': lambda x: x * config.base_count_multiplier,
                                           'scale': 'Blues', 'placement': 30})
        headers['gc content %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Spectral', 'placement': 40})
        headers['Q20 read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Reds', 'placement': 50})
        headers['Q30 base %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Oranges', 'placement': 60})
        headers['read length'].update({'suffix': 'bp', 'format': '{:,.0f}', 'scale': 'Greens', 'placement': 70})
        headers['no base call %'].update({'suffix': '%', 'format': '{:,.2f}', 'floor': 1, 'scale': 'BuGn'})
        # headers['bases sequenced (MB)'].update({'suffix': 'Mb', 'format': '{:,.2f}'})

        # general table: add read count and bases sequenced
        self.general_stats_addcols(self.qcdata,
                                   self.dict_ordered_subset(headers, ('read count', 'bases sequenced', 'gc content %')))

        log.info("Found {} reports".format(len(self.qcdata)))

        # write full data set to file
        self.write_data_file(self.qcdata, 'multiqc_readqc')

        # overview table with all values
        self.add_section(
            name='Overview',
            anchor='readqc-all',
            description='',
            plot=table.plot(self.qcdata,
                            headers,
                            pconfig={'namespace': 'ReadQC'})
        )
