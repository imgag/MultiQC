#!/usr/bin/env python

"""MultiQC module to parse output from VariantQC"""

from multiqc.modules.readqc import QcmlMultiqcModule
from multiqc.plots import table, bargraph

import logging

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(QcmlMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='VariantQC',
                                            anchor='variantqc',
                                            href="https://github.com/imgag/ngs-bits",
                                            info="calculates QC metrics based on variant lists.")

        # quality parameters from qcML with name, accession, description
        self.qcml = dict()
        # qc data for each sample
        self.qcdata = dict()
        # parse qcml files
        for f in self.find_log_files('variantqc', filecontents=True, filehandles=False):
            self.add_data_source(f)
            s_name = self.clean_s_name(f['s_name'], f['root'])
            self.qcdata[s_name] = self.parse_qcml(f['f'])

        # ignore samples if requested
        self.qcdata = self.ignore_samples(self.qcdata)

        # warn if no samples found
        if len(self.qcdata) == 0:
            raise UserWarning

        # prepare table headers, use name and description from qcML
        headers = {qp_key: {
            'namespace': "VariantQC",
            'title': qp_key,
            'description': qp_entry['description'],
        } for qp_key, qp_entry in self.qcml.items()}

        headers['variant count'].update({'format': '{:,.0f}', 'scale': 'Blues'})
        headers['known variants %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'YlGnBu'})
        headers['high-impact variants %'].update(
            {'suffix': '%', 'format': '{:,.2f}', 'min': 0, 'minRange': 10, 'ceiling': 10, 'scale': 'Reds'})
        headers['homozygous variants %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Purples'})
        headers['indel variants %'].update({'suffix': '%', 'format': '{:,.2f}', 'minRange': 20, 'ceiling': 20, 'scale': 'PuRd'})
        headers['transition/transversion ratio'].update({'format': '{:,.2f}', 'minRange': 5, 'ceiling': 5, 'scale': 'RdBu'})

        # general table: add read count and bases usable
        self.general_stats_addcols(self.qcdata,
                                   self.dict_ordered_subset(headers, (
                                       'variant count',
                                       'known variants %')))

        # write full data set to file
        self.write_data_file(self.qcdata, 'multiqc_variantqc')

        # table with general values
        self.add_section(
            name='Overview',
            anchor='variantqc-general',
            description='',
            plot=table.plot(self.qcdata,
                            self.dict_ordered_subset(headers, (
                                'variant count',
                                'high-impact variants %',
                                'homozygous variants %',
                                'indel variants %',
                                'transition/transversion ratio'
                            )),
                            pconfig={'namespace': 'VariantQC'})
        )

        # bar plot with variant count values
        self.add_section(
            name='Variant Count',
            anchor='variantqc-variant-count',
            description=self.make_description(['variant count']),
            plot=bargraph.plot(
                self.qcdata,
                self.dict_ordered_subset(headers, ('variant count',)),
                pconfig={
                    'namespace': 'VariantQC',
                    'id': 'variantqc-variant-count-plot',
                    'title': 'VariantQC: Variant Count',
                    'ylab': 'count',
                    'yDecimals': False,
                    'cpswitch': False,
                    'tt_decimals': 0,
                    'tt_suffix': '',
                    'tt_percentages': False
                }
            )
        )