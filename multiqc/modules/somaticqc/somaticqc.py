#!/usr/bin/env python

"""MultiQC module to parse output from SomaticQC"""

from multiqc.modules.readqc import QcmlMultiqcModule
from multiqc.plots import table, bargraph

import logging
import re

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(QcmlMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='SomaticQC',
                                            anchor='somaticqc',
                                            href="https://github.com/imgag/ngs-bits",
                                            info="calculates QC metrics based on tumor-normal pairs.")

        # quality parameters from qcML with name, accession, description
        self.qcml = dict()
        # qc data for each sample
        self.qcdata = dict()
        # parse qcml files
        for f in self.find_log_files('somaticqc', filecontents=True, filehandles=False):
            self.add_data_source(f)
            s_name = self.clean_s_name(f['s_name'], f['root'])
            # try to split Sample1-Sample2 names
            ms = re.match(r'([^-]+)-[^-]+', s_name)
            if (ms):
                s_name = ms.group(1)
            self.qcdata[s_name] = self.parse_qcml(f['f'])

        # ignore samples if requested
        self.qcdata = self.ignore_samples(self.qcdata)

        # warn if no samples found
        if len(self.qcdata) == 0:
            raise UserWarning

        # parse somatic variant rate value
        for s, kv in self.qcdata.items():
            try:
                kv['somatic variant rate'] = float(
                    re.sub(r'(low|moderate|high) \(([0-9.]+) var/Mb\)', '\2', kv['somatic variant rate']))
            except ValueError as e:
                kv.pop('somatic variant rate')

        # prepare table headers, use name and description from qcML
        headers = {qp_key: {
            'namespace': "SomaticQC",
            'title': qp_key,
            'description': qp_entry['description'],
        } for qp_key, qp_entry in self.qcml.items()}

        headers['sample correlation'].update({'format': '{:,.2f}', 'max': 1})
        headers['variant count'].update({'format': '{:,.0f}', 'title': 'variant count'})
        headers['somatic variant count'].update({'format': '{:,.0f}'})
        headers['known somatic variants %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100})
        headers['somatic indel %'].update({'suffix': '%', 'format': '{:,.2f}', 'minRange': 20, 'ceiling': 20})
        headers['somatic variant rate'].update(
            {'suffix': 'Variants/Mb', 'format': '{:,.2f}', 'min': 0, 'minRange': 10, 'ceiling': 10})

        try:
            headers['somatic transition/transversion ratio'].update({'format': '{:,.2f}', 'minRange': 5, 'ceiling': 5})
        except KeyError as e:
            pass

        try:
            headers['tumor content estimate'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100})
        except KeyError as e:
            pass

        # rename 'variant count' key to prevent duplicate ID with 'variant count' from VariantQC
        headers['variant count somaticqc'] = headers.pop('variant count')
        for s, kv in self.qcdata.items():
            kv['variant count somaticqc'] = kv.pop('variant count')

        # general table: add read count and bases usable
        self.general_stats_addcols(self.qcdata,
                                   self.dict_ordered_subset(headers, (
                                   'sample correlation', 'somatic variant count', 'known somatic variants %')))

        # write full data set to file
        self.write_data_file(self.qcdata, 'multiqc_somaticqc')

        # table with general values
        self.add_section(
            name='Overview',
            anchor='somaticqc-general',
            description='',
            plot=table.plot(self.qcdata,
                            self.dict_ordered_subset(headers, (
                                'sample correlation',
                                'variant count somaticqc',
                                'somatic variant count',
                                'known somatic variants %',
                                'somatic indel %',
                                'somatic transition/transversion ratio',
                                'somatic variant rate',
                                'tumor content estimate'
                            )),
                            pconfig={'namespace': 'SomaticQC'})
        )

        # bar plot with variant count values
        self.add_section(
            name='Somatic Variant Count',
            anchor='somaticqc-somatic-variant-count',
            description=self.make_description(['somatic variant count']),
            plot=bargraph.plot(
                self.qcdata,
                self.dict_ordered_subset(headers, ('somatic variant count',)),
                pconfig={
                    'namespace': 'SomaticQC',
                    'id': 'somaticqc-somatic-variant-count-plot',
                    'title': 'SomaticQC: Somatic Variant Count',
                    'ylab': 'count',
                    'yDecimals': False,
                    'cpswitch': False,
                    'tt_decimals': 0,
                    'tt_suffix': '',
                    'tt_percentages': False
                }
            )
        )