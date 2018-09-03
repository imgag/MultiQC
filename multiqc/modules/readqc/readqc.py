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

    def parse_qcml(self, qcml_contents):
        """Parse a qcML file and return key-value pairs from the quality parameter entries."""
        root = xml.etree.cElementTree.fromstring(qcml_contents)
        quality_parameters = dict()
        for qp in root.findall(".//{http://www.prime-xs.eu/ms/qcml}qualityParameter"):
            # skip n/a values
            if qp.attrib['value'].startswith('n/a'):
                continue

            # replace 'percentage' with '%'
            qp_name = re.sub(r' percentage$', ' %', qp.attrib['name'])

            try:
                quality_parameters[qp_name] = float(qp.attrib['value'])
            except ValueError:
                quality_parameters[qp_name] = qp.attrib['value']

            self.qcml[qp_name] = {'description': qp.attrib['description'],
                                  'accession': qp.attrib['accession']}
        return quality_parameters

    def make_description(self, keynames):
        """"Create description string from qcML quality parameter key name."""
        if len(keynames) == 1:
            desc = "{:s}".format(self.qcml[keynames[0]]['description'],
                                 self.qcml[keynames[0]]['accession'])
        else:
            desc = "<ul>" + "".join(
                ["<li>{:s}</li>".format(self.qcml[key]['description'],
                                        self.qcml[key]['accession']) for key in
                 keynames]) + "</ul>"
        return desc

    @staticmethod
    def dict_ordered_subset(d, ks):
        """Return subset of a dictionary as an OrderedDict object, ignoring non-existent keys."""
        od = collections.OrderedDict()
        for k in ks:
            try:
                od[k] = d[k]
            except KeyError as e:
                pass
        return od


class MultiqcModule(QcmlMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='ReadQC',
                                            anchor='readqc',
                                            href="https://github.com/imgag/ngs-bits",
                                            info="calculates QC metrics on unprocessed NGS reads.")

        # quality parameters from qcML with name, accession, description
        self.qcml = dict()
        # qc data for each sample
        self.qcdata = dict()
        # parse qcml files
        for f in self.find_log_files('readqc', filecontents=True, filehandles=False):
            self.add_data_source(f)
            s_name = self.clean_s_name(f['s_name'], f['root'])
            self.qcdata[s_name] = self.parse_qcml(f['f'])

        # ignore samples if requested
        self.qcdata = self.ignore_samples(self.qcdata)

        # warn if no samples found
        if len(self.qcdata) == 0:
            raise UserWarning

        # add bases sequenced key, derived from bases sequenced (MB)
        self.qcml.pop('bases sequenced (MB)')
        self.qcml['bases sequenced'] = dict()
        self.qcml['bases sequenced']['description'] = 'Bases sequenced in total.'
        for s, kv in self.qcdata.items():
            kv['bases sequenced'] = kv['bases sequenced (MB)'] * 1e6
            kv.pop('bases sequenced (MB)')

        # prepare table headers, use name and description from qcML
        headers = {qp_key: {
            'namespace': "ReadQC",
            'title': qp_key,
            'description': qp_entry['description'],
        } for qp_key, qp_entry in self.qcml.items()}
        headers['Q20 read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100})
        headers['Q30 base %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100})
        headers['gc content %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100})
        headers['no base call %'].update({'suffix': '%', 'format': '{:,.2f}', 'floor': 1})
        # headers['bases sequenced (MB)'].update({'suffix': 'Mb', 'format': '{:,.2f}'})
        headers['bases sequenced'].update({'suffix': config.base_count_prefix, 'format': '{:,.2f}',
                                           'modify': lambda x: x * config.base_count_multiplier})
        headers['read count'].update({'suffix': config.read_count_prefix, 'format': '{:,.2f}',
                                      'modify': lambda x: x * config.read_count_multiplier})
        headers['read length'].update({'suffix': 'bp', 'format': '{:,.0f}'})

        # general table: add read count and bases sequenced
        self.general_stats_addcols(self.qcdata,
                                   self.dict_ordered_subset(headers, ('read count', 'bases sequenced', 'gc content %')))

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
