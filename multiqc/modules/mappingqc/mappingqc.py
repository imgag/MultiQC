#!/usr/bin/env python

"""MultiQC module to parse output from MappingQC"""

from multiqc.modules.readqc import QcmlMultiqcModule
from multiqc import config
from multiqc.plots import table, bargraph

import logging

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(QcmlMultiqcModule):

    def __init__(self):
        super(MultiqcModule, self).__init__(name='MappingQC',
                                            anchor='mappingqc',
                                            href="https://github.com/imgag/ngs-bits",
                                            info="calculates QC metrics based on mapped NGS reads.")

        # quality parameters from qcML with name, accession, description
        self.qcml = dict()
        # qc data for each sample
        self.qcdata = dict()
        # parse qcml files
        for f in self.find_log_files('mappingqc', filecontents=True, filehandles=False):
            self.add_data_source(f)
            s_name = self.clean_s_name(f['s_name'], f['root'])
            self.qcdata[s_name] = self.parse_qcml(f['f'])

        # ignore samples if requested
        self.qcdata = self.ignore_samples(self.qcdata)

        # warn if no samples found
        if len(self.qcdata) == 0:
            raise UserWarning

        # add bases usable key, derived from bases usable (MB)
        self.qcml.pop('bases usable (MB)')
        self.qcml['bases usable'] = dict()
        self.qcml['bases usable']['description'] = 'Bases sequenced in total.'
        for s, kv in self.qcdata.items():
            kv['bases usable'] = kv['bases usable (MB)'] * 1e6
            kv.pop('bases usable (MB)')

        # prepare table headers, use name and description from qcML
        headers = {qp_key: {
            'namespace': "MappingQC",
            'title': qp_key,
            'description': qp_entry['description'],
        } for qp_key, qp_entry in self.qcml.items()}

        headers['trimmed base %'].update({'suffix': '%', 'format': '{:,.2f}', 'floor': 1, 'scale': 'PuBu'})
        headers['clipped base %'].update({'suffix': '%', 'format': '{:,.2f}', 'floor': 1, 'scale': 'PuRd'})
        headers['mapped read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Reds'})
        headers['bases usable'].update({'suffix': config.base_count_prefix, 'format': '{:,.2f}',
                                        'modify': lambda x: x * config.base_count_multiplier,
                                        'scale': 'Greens'})
        # always available, even without target file
        headers['on-target read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'Purples'})

        # only available if duplicates marked
        try:
            headers['duplicate read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'YlOrRd'})
        except KeyError:
            pass

        # only available if paired-end
        try:
            headers['properly-paired read %'].update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'GnBu'})
            headers['insert size'].update({'suffix': 'bp', 'format': '{:,.2f}', 'scale': 'RdYlGn'})
        except KeyError:
            pass

        # only available if human
        try:
            headers['SNV allele frequency deviation'].update(
                {'suffix': '', 'format': '{:,.2f}', 'floor': 0, 'ceiling': 10, 'minRange': 10, 'scale': 'Greys'})
        except KeyError:
            pass

        # only available if target file provided
        coverage_values = (10, 20, 30, 50, 100, 200, 500)
        try:
            headers['target region read depth'].update({'suffix': 'x', 'format': '{:,.2f}'})
            for x in coverage_values:
                headers['target region {:d}x %'.format(x)]. \
                    update({'suffix': '%', 'format': '{:,.2f}', 'max': 100, 'scale': 'YlGn'})
        except KeyError:
            pass

        # general table: add read count and bases usable
        self.general_stats_addcols(self.qcdata,
                                   self.dict_ordered_subset(headers, (
                                       'bases usable',
                                       'mapped read %',
                                       'on-target read %',
                                       'target region read depth')))

        # write full data set to file
        self.write_data_file(self.qcdata, 'multiqc_mappingqc')

        # table with general values
        self.add_section(
            name='Overview',
            anchor='mappingqc-general',
            description='',
            plot=table.plot(self.qcdata,
                            self.dict_ordered_subset(headers, (
                                'bases usable',
                                'on-target read %',
                                'mapped read %',
                                'properly-paired read %',
                                'trimmed base %',
                                'clipped base %',
                                'duplicate read %',
                                'insert size',
                                'SNV allele frequency deviation'
                            )),
                            pconfig={'namespace': 'MappingQC'})
        )

        if 'target region 10x %' in headers.keys():
            # table with coverage values
            self.add_section(
                name='Coverage',
                anchor='mappingqc-coverage',
                description='',
                plot=table.plot(self.qcdata,
                                self.dict_ordered_subset(headers, (
                                    'target region read depth',
                                    'target region 10x %',
                                    'target region 20x %',
                                    'target region 30x %',
                                    'target region 50x %',
                                    'target region 100x %',
                                    'target region 200x %',
                                    'target region 500x %'
                                )),
                                pconfig={'namespace': 'MappingQC'})
            )

            # bar plot with sequencing depth values
            self.add_section(
                name='Sequencing Depth',
                anchor='mappingqc-read-depth',
                description=self.make_description(['target region read depth']),
                plot=bargraph.plot(
                    self.qcdata,
                    self.dict_ordered_subset(headers, ('target region read depth',)),
                    pconfig={
                        'namespace': 'MappingQC',
                        'id': 'mappingqc-read-depth-plot',
                        'title': 'MappingQC: Target Region Sequencing Depth',
                        'ylab': 'coverage',
                        'cpswitch': False,
                        'tt_decimals': 2,
                        'tt_suffix': 'x',
                        'tt_percentages': False
                    }
                )
            )

            # bar plot with coverage values
            self.add_section(
                name='Target Coverage',
                anchor='mappingqc-target-coverage',
                description='',
                plot=bargraph.plot(
                    [self.qcdata] * len(coverage_values),
                    [{s: headers[s]} for s in ['target region {:d}x %'.format(x) for x in coverage_values]],
                    pconfig={
                        'namespace': 'MappingQC',
                        'id': 'mappingqc-target-coverage-plot',
                        'title': 'MappingQC: Target Coverage Percentage',
                        'ylab': 'target coverage percentage',
                        'cpswitch': False,
                        'data_labels': ['{:d}x coverage %'.format(x) for x in coverage_values],
                        'ymin': 0,
                        'ymax': 100,
                        'use_legend': False,
                        'tt_decimals': 2,
                        'tt_suffix': '%',
                        'tt_percentages': False
                    }
                )
            )
