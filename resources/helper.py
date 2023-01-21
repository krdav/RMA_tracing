import pandas as pd
import numpy as np
import copy as cp
import warnings


class PeakData:
    def __init__(self, name, sample_info, params):
        self.name = name
        self.hydrogen_mass = 1.007825031898
        self.electron_mass = 5.485799e-4
        self.sample_info = cp.deepcopy(sample_info)
        self.params = cp.deepcopy(params)
        # Dummy for peak data:
        self.peak_data_pos = None
        self.peak_data_neg = None
        # Dummy for column names containing peak area:
        self.area_colnames_pos = None
        self.area_colnames_neg = None

        # For each label make a dictionary entry with all the
        # tables related to the peak pair from that label:
        self.labels = {self.sample_info[sample]['label'] for sample in self.sample_info if self.sample_info[sample]['type'] != 'blank' and self.sample_info[sample]['label'] != 'None'}
        self.sample_label = dict()
        self.labelled_names = set()
        for label in self.labels:
            self.sample_label[label] = [self.sample_info[sample]['name'] for sample in self.sample_info if self.sample_info[sample]['label'] == label and self.sample_info[sample]['type'] != 'blank']
            self.labelled_names = self.labelled_names.union(self.sample_label[label])

        self.label_pairs = dict()
        # Iterate over labels:
        for label in sorted(self.labels):
            self.label_pairs[label] = dict()
            self.label_pairs[label]['label_colnames'] = self.sample_label[label]
            for polarity in ['pos', 'neg']:
                self.label_pairs[label][polarity] = dict()
                self.label_pairs[label][polarity]['peak_pair_area_parent'] = None # Parent peak area
                self.label_pairs[label][polarity]['peak_pair_area_heavy'] = None # Labelled peak area
                self.label_pairs[label][polarity]['peak_pair_labelp'] = None # Labelling percent
                self.label_pairs[label][polarity]['area_ratio_mask'] = None # Filter based on self.params e.g. mass different, retention time difference etc.
                self.label_pairs[label][polarity]['peak_pair_corr'] = None # Matrix with accross sample correlation coefficients between peak areas

    # Read peak data according to polarity:
    def read_peaks(self, datafile, polarity):
        if polarity == 'pos':
            self.peak_data_pos, self.area_colnames_pos = self.__peak_reader(datafile, polarity)
        elif polarity == 'neg':
            self.peak_data_neg, self.area_colnames_neg = self.__peak_reader(datafile, polarity)
        else:
            raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    # Internal function to read peaks in either polarity:
    def __peak_reader(self, datafile, polarity):
        csv = pd.read_excel(datafile)
        csv.columns = [n.lower().replace(' ', '_') for n in csv.columns]
        # Assert that weight and retention time can be used as unique identifier:
        wr_list = [(w, r) for w, r, in zip(csv['molecular_weight'], csv['rt_[min]'])]
        assert(len(wr_list) == len(set(wr_list)))
        # Insert (MW, RT) as peak ID into the dataframe:
        csv.insert(0, 'MW', csv['molecular_weight'])
        csv.insert(0, 'RT', csv['rt_[min]'])
        csv.insert(0, 'peak_id', pd.Series(wr_list, index=csv.index))
        # Extract the columns with area information:
        col_sele = ['peak_id', 'RT', 'MW', 'name']
        col_sele.extend([col for col in csv.columns if 'area' in col])
        peak_data = csv.loc[:, col_sele]

        # Get sample names from sample description:
        colname_list = list()
        area_colnames = list()
        area_colnames_type = list()
        area_colnames_order = list()
        for col in peak_data.columns:
            if 'area:' in col:
                colname = False
                for sample in self.sample_info.keys():
                    if colname is False and sample.lower() in col.lower():
                        colname = self.sample_info[sample]['name']
                        order = self.sample_info[sample]['order']
                        sample_type = self.sample_info[sample]['type']
                    elif colname is True and sample.lower() in col.lower():
                        raise Exception('Column "{}" maps to multiple names in sample description. Please correct.'.format(col))
                if colname is False:
                        raise Exception('Column name "{}" not found. Check sample description.'.format(col))
                colname_list.append(colname)
                area_colnames.append(colname)
                area_colnames_type.append(sample_type)
                area_colnames_order.append(order)
            else:
                colname_list.append(col)
        peak_data.columns = colname_list

        # Get the area column names in order:
        sample_area_colnames = [(nam, order) for nam, smtyp, order in zip(area_colnames, area_colnames_type, area_colnames_order) if smtyp == 'sample']
        sample_area_colnames = [val[0] for val in sorted(sample_area_colnames, key=lambda x: x[1])]
        blank_area_colnames = [(nam, order) for nam, smtyp, order in zip(area_colnames, area_colnames_type, area_colnames_order) if smtyp == 'blank']
        blank_area_colnames = [val[0] for val in sorted(blank_area_colnames, key=lambda x: x[1])]

        # Filter peak dataframe:
        print('Running peak filtering for polarity: {}'.format(polarity))
        peak_data = self.__run_all_peak_filters(peak_data, sample_area_colnames, blank_area_colnames)

        # Drop blank peak areas, then sort peak areas columns
        # for samples using join:
        peak_data = peak_data.drop(labels=blank_area_colnames, axis=1)
        peak_data = peak_data.loc[:, ~peak_data.columns.isin(sample_area_colnames)].join(peak_data.loc[:, sample_area_colnames])
        # Sort according to molecular weight then return:
        peak_data = peak_data.sort_values(by='MW', ascending=True)
        peak_data.reset_index(drop=True, inplace=True)

        return(peak_data, sample_area_colnames)

    # Apply all filters on the peak data:
    def __run_all_peak_filters(self, peak_data, sample_area_colnames, blank_area_colnames):
        '''
        Function to run sequential filtering on peaks.
        '''
        count_before = len(peak_data)
        # Min area filter:
        peak_data = self.min_area_peak_filter(peak_data, sample_area_colnames)
        N_min_area = count_before - len(peak_data)
        # Min MW filter:
        peak_data = self.min_MW_peak_filter(peak_data, self.params['min_MW'])
        N_min_MW = count_before - len(peak_data) - N_min_area
        # Blank filter:
        peak_data = self.blank_peak_filter(peak_data, sample_area_colnames, blank_area_colnames, self.params['min_fold_blank'])
        N_blank = count_before - len(peak_data) - N_min_area - N_min_MW
        # Merge closely related:
        # Apply RT and correlation criterium separately.
        peak_data = self.merge_peak_filter(peak_data, sample_area_colnames, self.params['merge_ppm_tol'], self.params['merge_RT_tol'], -1)
        peak_data = self.merge_peak_filter(peak_data, sample_area_colnames, self.params['merge_ppm_tol'], self.params['merge_RT_tol']*2, self.params['merge_corr_tol'])
        N_merge = count_before - len(peak_data) - N_min_area - N_min_MW - N_blank
        count_after = len(peak_data)
        print('Filtered {} peaks out based on.\n\
Minimum peak area: {}\n\
Minimum molecular weight: {}\n\
Minimum fold over blank: {}\n\
Merged closely related peaks: {}\n\
{} peaks left.\n'.format(count_before - count_after, N_min_area, N_min_MW, N_blank, N_merge, count_after))
        return(peak_data)

    # Filter peaks with area smaller than the cutoff:
    def min_area_peak_filter(self, peak_data, area_colnames):
        '''
        Enforce a minimum peak area.
        For a given peak the peak area at least one sample has to be above this value.
        '''
        # Global minimum area:
        area_max = peak_data.loc[:, area_colnames].max(axis=1)
        min_mask = (area_max > self.params['min_area'])
        # Mininum area for any of the labelled samples:
        labelled_columns = [cn for cn in area_colnames if cn in self.labelled_names]
        area_max_label = peak_data.loc[:, labelled_columns].max(axis=1)
        min_mask_label = (area_max_label > self.params['min_area_label'])

        return(peak_data[min_mask&min_mask_label])

    # Filter peaks with molecular weight smaller than the cutoff:
    def min_MW_peak_filter(self, peak_data, min_MW):
        mask = peak_data['MW'] > min_MW
        return(peak_data[mask])

    # Filter peaks with maximum peak area less than X fold higher
    # than max peak area for blank samples:
    def blank_peak_filter(self, peak_data, sample_area_colnames, blank_area_colnames, min_fold):
        sample_area_max = peak_data.loc[:, sample_area_colnames].max(axis=1)
        blank_area_max = peak_data.loc[:, blank_area_colnames].max(axis=1)
        mask = (sample_area_max / blank_area_max) > min_fold
        return(peak_data[mask])

    # Merge peaks with small deviation in m/z.
    def merge_peak_filter(self, peak_data, area_colnames, ppm_tol, RT_tol, corr_tol):
        '''
        This filter is merging peaks that are likely from the same compound.
        It works by searching for peaks that fulfill two criteria:
        1) Being within a maximum mass difference (defined by a ppm value).
        2) Being within a maximum retention time difference,
        OR Being within a maximum retention time difference x2 AND
        having a minimum correlation coefficient between the two peak areas.

        Peaks are merged by taken the sum of the peak areas and keeping
        the molecular mass and retention time from the peak with the 
        largest sum of peak areas.
        '''
        MW = peak_data['MW'].values
        RT = peak_data['RT'].values
        area = peak_data.loc[:, area_colnames].values.astype(float)

        # Find peaks that fulfill the two criteria:
        merge_col = list()
        for i in range(len(peak_data)):
            # Retention time criterium:
            RT_diff_mask = np.abs(RT[i] - RT) <= RT_tol
            # Mass shift criterium:
            MW_tol = MW[i] * 1e-6 * ppm_tol
            MW_diff_mask = np.abs(MW[i] - MW) <= MW_tol
            # Peak area correlation criterium:
            corr_mask = vcorrcoef(area, area[i]) > corr_tol

            # Make a mask and store peak pairs:
            mask = MW_diff_mask & RT_diff_mask & corr_mask
            idx_tup = tuple(sorted(np.where(mask)[0]))
            merge_col.append(idx_tup)

        # Add a "merge" column and split the dataframe
        # into peak areas and other data: 
        peak_data['merge'] = merge_col
        peak_data_area = peak_data.loc[:, ['merge'] + area_colnames].copy()
        peak_data_rest = peak_data.loc[:, ~peak_data.columns.isin(area_colnames)].copy()
        # Add the sum of area for each peak, as a way of sorting:
        peak_data_rest['area_sum'] = peak_data_area.sum(axis=1, numeric_only=True).values

        # Group dataframes by the merge column.
        # Take the sum of the areas to the merged:
        peak_data_area_grouped = peak_data_area.groupby('merge').sum().reset_index()
        # Extract the molecular weight and retention time
        # of the peak with the largest sum of peak areas:
        largest_area_sum_mask = peak_data_rest.groupby(['merge'])['area_sum'].transform(max) == peak_data_rest['area_sum']
        peak_data_rest_grouped = peak_data_rest[largest_area_sum_mask].copy()

        # Merge back the area dataframe with the rest
        # and remove ancillary columns:
        peak_data_merged = peak_data_rest_grouped.merge(peak_data_area_grouped)
        peak_data_merged = peak_data_merged.drop(labels=['merge', 'area_sum'], axis=1)

        return(peak_data_merged)

    def annotate_known_peaks(self, known_fnam, formula2mass):
        '''
        Search for and add compound name to peaks
        that have a matching known compound in
        the input list.
        '''
        df_known = pd.read_csv(known_fnam, sep='\t')
        # Convert the "Formula" column to an exact mass:
        df_known['Mass'] = [formula2mass(f) for f in df_known['Formula']]

        # Annotate in both polarities:
        for polarity in ['pos', 'neg']:
            # Skip if no data for this polarity,
            # otherwise add annotations:
            if polarity == 'pos':
                if self.area_colnames_pos == None:
                    continue
                self.peak_data_pos = self.__annotate_known_polarity(df_known, self.peak_data_pos)
            elif polarity == 'neg':
                if self.area_colnames_neg == None:
                    continue
                self.peak_data_neg = self.__annotate_known_polarity(df_known, self.peak_data_neg)

    def __annotate_known_polarity(self, df_known, peak_data):
        '''Add annotations in a specified polarity.'''

        # Store annotations in this dictionary:
        anno_dict = dict()

        # Search for each known compound in the peak data:
        for index, row in df_known.iterrows():
            mass_error = row['MW_ppm_tol'] * row['Mass'] * 1e-6
            RT_mask = ((row['RT'] + row['RT_tol']) > peak_data['RT']) & ((row['RT'] - row['RT_tol']) < peak_data['RT'])
            MW_mask = ((row['Mass'] + mass_error) > peak_data['MW']) & ((row['Mass'] - mass_error) < peak_data['MW'])
            mask = RT_mask & MW_mask
            # Store annotations in index based dictionary:
            if mask.sum() > 0:
                for j in np.where(mask)[0]:
                    # If multiple known compounds are found,
                    # store the one with closest retention time:
                    if j in anno_dict:
                        old_RT_diff = np.abs(anno_dict[j][1] - peak_data['RT'][j])
                        new_RT_diff = np.abs(row['RT'] - peak_data['RT'][j])
                        if old_RT_diff > new_RT_diff:
                            anno_dict[j] = (row['Name'], row['RT'])
                    else:
                        anno_dict[j] = (row['Name'], row['RT'])

        # Convert annotation dictionary to list,
        # and add to peak dataframe:
        anno_list = [anno_dict[i][0] if i in anno_dict else None for i in range(len(peak_data))]
        peak_data['known_anno'] = anno_list

        return(peak_data)

    def remove_blacklist_peaks(self, blacklist, polarity='both'):
        ''' Remove blacklisted peaks from the peak data.'''
        if polarity in ['pos', 'both']:
            count_before = len(self.peak_data_pos)
            MW = self.peak_data_pos.loc[:, 'MW'].values
            RT = self.peak_data_pos.loc[:, 'RT'].values
            blacklist_mask = MW > 0 # dummy mask, all True
            for peak in blacklist['pos']:
                MW_i = blacklist['pos'][peak]['MW']
                RT_i = blacklist['pos'][peak]['RT']
                # Retention time criterium:
                RT_diff_mask = np.abs(RT_i - RT) <= blacklist['pos'][peak]['RT_tol']
                # Mass shift criterium:
                MW_tol = MW_i * 1e-6 * blacklist['pos'][peak]['MW_ppm_tol']
                MW_diff_mask = np.abs(MW_i - MW) <= MW_tol
                blacklist_mask = blacklist_mask & ~(RT_diff_mask & MW_diff_mask)
            self.peak_data_pos = self.peak_data_pos[blacklist_mask].reset_index(drop=True, inplace=False)
            count_after = len(self.peak_data_pos)
            print('Blacklist filter in positive polarity filtered {} peaks out. {} peaks left.'.format(count_before - count_after, count_after))

        if polarity in ['neg', 'both']:
            count_before = len(self.peak_data_neg)
            MW = self.peak_data_neg.loc[:, 'MW'].values
            RT = self.peak_data_neg.loc[:, 'RT'].values
            blacklist_mask = MW > 0 # dummy mask, all True
            for peak in blacklist['neg']:
                MW_i = blacklist['neg'][peak]['MW']
                RT_i = blacklist['neg'][peak]['RT']
                # Retention time criterium:
                RT_diff_mask = np.abs(RT_i - RT) <= blacklist['neg'][peak]['RT_tol']
                # Mass shift criterium:
                MW_tol = MW_i * 1e-6 * blacklist['neg'][peak]['MW_ppm_tol']
                MW_diff_mask = np.abs(MW_i - MW) <= MW_tol
                blacklist_mask = blacklist_mask & ~(RT_diff_mask & MW_diff_mask)
            self.peak_data_neg = self.peak_data_neg[blacklist_mask].reset_index(drop=True, inplace=False)
            count_after = len(self.peak_data_neg)
            print('Blacklist filter in negative polarity filtered {} peaks out. {} peaks left.'.format(count_before - count_after, count_after))
        if polarity not in ['pos', 'neg', 'both']:
            raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    def pick_ratio(self, known_fnam, polarity, formula2mass, label, label_str):
        '''
        Search for known peak pairs and calculate statistics
        to help adjust filtering parameters.
        The input known labelled compounds is provided as a tab
        separated file.
        The polarity for the search is specified.
        A function ("formula2mass") is specified from an Isotope
        object to convert formula to mass.
        The label in specified.
        The "label string" is specified as a list of two strings,
        indicating the nominal mass shift for each isotope of the compound.
        Example: ['m', m+4] to compare the parent unlabelled compound
        to the nominally 4 Da mass shifted labelled compound.
        '''

        try:
            labeled_columns = self.sample_label[label]
        except KeyError:
            raise Exception('Label not defined: {}'.format(label))

        if polarity == 'pos':
            peak_data = self.peak_data_pos
            area_colnames = self.area_colnames_pos
        elif polarity == 'neg':
            peak_data = self.peak_data_neg
            area_colnames = self.area_colnames_neg
        else:
            raise Exception('Polarity not valid: {}'.format(polarity))


        df_known = pd.read_csv(known_fnam, sep='\t')
        # Convert the "Formula" column to an exact mass:
        df_known['Mass'] = [formula2mass(f) for f in df_known['Formula']]
        known_compounds = sorted(list(set(df_known['Name'].values)))
        known_idx = [n + ' ({})'.format(f) for n, f in zip(df_known['Name'], df_known['Label'])]
        perc_lab = pd.DataFrame(columns=area_colnames, index=known_idx)

        # For each compound:
        for compound in known_compounds:
            dp_df = df_known[df_known['Name'] == compound]
            area = {label_str[0]:[0], label_str[1]:[0]}
            # Search in the peaks for each labelled compound
            # i.e. for each isotope: 
            for index, row in dp_df.iterrows():
                mass_error = row['MW_ppm_tol'] * row['Mass'] * 1e-6
                RT_mask = ((row['RT'] + row['RT_tol']) > peak_data['RT']) & ((row['RT'] - row['RT_tol']) < peak_data['RT'])
                MW_mask = ((row['Mass'] + mass_error) > peak_data['MW']) & ((row['Mass'] - mass_error) < peak_data['MW'])
                mask = RT_mask & MW_mask
                # The search might yield multiple peaks,
                # in such cases only use the peak with the largest peak area:
                mask_area_max = peak_data[area_colnames].sum(1) == peak_data[mask][area_colnames].sum(1).max()
                mask = mask & mask_area_max
                if sum(mask) == 1:
                    area[row['Label']] = peak_data[mask][area_colnames].values[0]

            # Find the percent labelling:
            if sum(area[label_str[0]]) > 0 and sum(area[label_str[1]]) > 0:
                area_sum = area[label_str[0]] + area[label_str[1]]
                perc_m = area[label_str[0]] / area_sum
                perc_mX = area[label_str[1]] / area_sum

                idx = compound + ' ({})'.format(label_str[0])
                perc_lab.loc[idx] = perc_m

                idx = compound + ' ({})'.format(label_str[1])
                perc_lab.loc[idx] = perc_mX

        x_labeled = perc_lab.loc[:, perc_lab.columns.isin((labeled_columns))].dropna().astype('float64')
        desc = x_labeled.T.describe(include='all')

        return(desc)

    # Find peak pairs according to polarity:
    def find_pairs(self, polarity):
        for label in sorted(self.labels):            
            if polarity == 'pos':
                self.label_pairs[label]['pos']['peak_pair_area_parent'], \
                self.label_pairs[label]['pos']['peak_pair_area_heavy'], \
                self.label_pairs[label]['pos']['peak_pair_labelp'], \
                self.label_pairs[label]['pos']['area_ratio_mask'], \
                self.label_pairs[label]['pos']['peak_pair_corr'] = \
                self.__find_peak_pairs_per_label(self.peak_data_pos, \
                                                 self.area_colnames_pos, \
                                                 polarity, label)
            elif polarity == 'neg':
                self.label_pairs[label]['neg']['peak_pair_area_parent'], \
                self.label_pairs[label]['neg']['peak_pair_area_heavy'], \
                self.label_pairs[label]['neg']['peak_pair_labelp'], \
                self.label_pairs[label]['neg']['area_ratio_mask'], \
                self.label_pairs[label]['neg']['peak_pair_corr'] = \
                self.__find_peak_pairs_per_label(self.peak_data_neg, \
                                                 self.area_colnames_neg, \
                                                 polarity, label)
            else:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    # Internal function to find peak pairs in either pos/neg polarity:
    def __find_peak_pairs_per_label(self, peak_data, area_colnames, polarity, label):
        '''
        In this function peak pairs are searched for by making an all-against-all comparison.
        Pairs have to fulfill the requirements of a maximum difference in retention time and
        a minimum proximity to the theoretical mass shift by incorporating the label.
        After this, the peak pairs are further filtered based on percentage label incorporation
        and can be filtered in many other ways.
        '''
        # Columns in peak pair table:
        pair_info_mask = ['pair_id', 'MW_parent', 'RT_parent', 'MW_heavy', 'RT_heavy',
                          'polarity', 'label', 'name', 'RT_diff', 'MW_ppm_diff', 'known_anno']
        pair_columns = cp.deepcopy(pair_info_mask)
        # Add area to the peak pair table:
        pair_columns.extend(area_colnames)
        # Make a table for the area of parent/heavy compound:
        peak_pair_area_parent = pd.DataFrame(columns = pair_columns)
        peak_pair_area_heavy = pd.DataFrame(columns = pair_columns)

        peak_pair_set = set()
        ppm_cutoff = self.params['pair_ppm_tol']
        RT_tol = self.params['pair_RT_tol']
        MW = peak_data.loc[:, 'MW'].values
        RT = peak_data.loc[:, 'RT'].values
        peak_id = peak_data.loc[:, 'peak_id'].values
        MW_shift = self.params['MW_shift'][label] # The mass shift is defined by the label
        # Iterate over peaks to match on MW/RT criteria:
        for i in range(len(peak_data)):
            # Retention time criterium:
            RT_diff_mask = np.abs(RT[i] - RT) <= RT_tol
            # Mass shift criterium:
            MW_tol = MW[i] * 1e-6 * ppm_cutoff
            MW_diff_high = MW_shift + MW_tol
            MW_diff_low = MW_shift - MW_tol
            MW_diff_higer_mask = ((MW - MW[i]) <= MW_diff_high) & ((MW - MW[i]) >= MW_diff_low)

            # Make a mask and store peak pairs:
            mask = (RT_diff_mask & MW_diff_higer_mask)
            if mask.sum() > 0:
                for j in np.where(mask)[0]:
                    # Because df is sorted by MW,
                    # i is parent, j is heavy
                    CD_name = peak_data.loc[i, 'name']
                    known_anno = peak_data.loc[i, 'known_anno']
                    pair_id = (peak_id[i], peak_id[j], polarity, label)
                    assert(pair_id not in peak_pair_set)
                    peak_pair_set.add(pair_id)

                    # Add data to each peak pair table:
                    tmp_df = pd.DataFrame(columns = pair_columns)
                    RT_diff = np.abs(RT[i] - RT[j])
                    MW_ppm_diff = np.abs((MW[i] + MW_shift) - MW[j]) / MW[i] * 1e6
                    tmp_df.loc[i, pair_info_mask] = np.array([pair_id, MW[i], RT[i], MW[j], RT[j], polarity, label, CD_name, RT_diff, MW_ppm_diff, known_anno], dtype=object)
                    tmp_df.loc[i, area_colnames] = peak_data.loc[i, area_colnames]
                    peak_pair_area_parent = pd.concat((peak_pair_area_parent, tmp_df), ignore_index=True)
                    tmp_df = pd.DataFrame(columns = pair_columns)
                    tmp_df.loc[j, area_colnames] = peak_data.loc[j, area_colnames]
                    tmp_df.loc[j, pair_info_mask] = np.array([pair_id, MW[i], RT[i], MW[j], RT[j], polarity, label, CD_name, RT_diff, MW_ppm_diff, known_anno], dtype=object)
                    peak_pair_area_heavy = pd.concat((peak_pair_area_heavy, tmp_df), ignore_index=True)

        # Reset index for peak pairs:
        peak_pair_area_parent.reset_index(drop=True, inplace=True)
        peak_pair_area_heavy.reset_index(drop=True, inplace=True)

        # Calculate the percent labelled for the peak pairs:
        peak_pair_labelp = peak_pair_area_heavy.loc[:, area_colnames] / (peak_pair_area_heavy.loc[:, area_colnames] + peak_pair_area_parent.loc[:, area_colnames])
        
        # Make a mask of the peak area ratio criteria, including min_area:
        area_ratio_cutoff = self.params['area_ratio_cutoff'][label][0]
        area_ratio_mask = ((peak_pair_labelp >= area_ratio_cutoff[0]) &  # ratio must be above min cutoff
                          (peak_pair_labelp <= area_ratio_cutoff[1]) &  # ratio must be below max cutoff
                          (peak_pair_area_parent.loc[:, area_colnames] > self.params['min_area'])) # peak area must be above the minimum
        # If multiple peak area ratio criteria, add the rest here:
        for area_ratio_cutoff in self.params['area_ratio_cutoff'][label][1:]:
            area_ratio_mask = (area_ratio_mask |  # already in the mask
                              ((peak_pair_labelp >= area_ratio_cutoff[0]) &  # ratio must be above min cutoff
                               (peak_pair_labelp <= area_ratio_cutoff[1])))  # ratio must be below max cutoff

        # Add pair info:
        pair_info_df = peak_pair_area_heavy.loc[:, pair_info_mask]
        peak_pair_labelp = pd.concat([pair_info_df, peak_pair_labelp], axis=1, sort=False)
        area_ratio_mask = pd.concat([pair_info_df, area_ratio_mask], axis=1, sort=False)

        # Find and drop rows with insufficient number of samples passing the area ratio criterium:
        area_ratio_drop = list()
        for i in range(len(peak_pair_area_parent)):
            cols_with_label = area_ratio_mask.columns.isin((self.label_pairs[label]['label_colnames']))
            if sum(area_ratio_mask.loc[i, cols_with_label]) < self.params['pair_min_area']:
                area_ratio_drop.append(i)

        peak_pair_area_parent = peak_pair_area_parent.drop(labels=area_ratio_drop, axis=0).reset_index(drop=True, inplace=False)
        peak_pair_area_heavy = peak_pair_area_heavy.drop(labels=area_ratio_drop, axis=0).reset_index(drop=True, inplace=False)
        peak_pair_labelp = peak_pair_labelp.drop(labels=area_ratio_drop, axis=0).reset_index(drop=True, inplace=False)
        area_ratio_mask = area_ratio_mask.drop(labels=area_ratio_drop, axis=0).reset_index(drop=True, inplace=False)

        # Make dataframe with peak pair area correlations:
        a = peak_pair_area_parent.sort_values(by='RT_parent', ascending=True).loc[:, area_colnames]
        a = a.astype(float)
        peak_pair_corr = a.T.corr().abs()
        pair_info_df = peak_pair_area_parent.sort_values(by='RT_parent', ascending=True).loc[:, pair_info_mask]
        peak_pair_corr = pd.concat([pair_info_df, peak_pair_corr], axis=1, sort=False)

        return(peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr)


    # Find adducts for peak pairs according to polarity:
    def flag_adducts(self, adducts_fnam, polarity='both'):   
        for label in sorted(self.labels):        
            if polarity in ['pos', 'both']:
                # Find adducts:
                adduct_flags = self.__flag_adducts_per_label(self.label_pairs[label]['pos']['peak_pair_area_parent'], self.area_colnames_pos, 'pos', adducts_fnam)
                corr_sort_order_rv = self.label_pairs[label]['pos']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                adduct_flags_corr_sorted = [t[0] for t in sorted(zip(adduct_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['pos']['peak_pair_area_parent'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['pos']['peak_pair_area_heavy'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['pos']['peak_pair_labelp'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['pos']['area_ratio_mask'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['pos']['peak_pair_corr'].insert(10, "Adducts", adduct_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['pos']['peak_pair_area_parent']["Adducts"] = adduct_flags
                    self.label_pairs[label]['pos']['peak_pair_area_heavy']["Adducts"] = adduct_flags
                    self.label_pairs[label]['pos']['peak_pair_labelp']["Adducts"] = adduct_flags
                    self.label_pairs[label]['pos']['area_ratio_mask']["Adducts"] = adduct_flags
                    self.label_pairs[label]['pos']['peak_pair_corr']["Adducts"] = adduct_flags_corr_sorted

            if polarity in ['neg', 'both']:
                # Find adducts:
                adduct_flags = self.__flag_adducts_per_label(self.label_pairs[label]['neg']['peak_pair_area_parent'], self.area_colnames_neg, 'neg', adducts_fnam)
                corr_sort_order_rv = self.label_pairs[label]['neg']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                adduct_flags_corr_sorted = [t[0] for t in sorted(zip(adduct_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['neg']['peak_pair_area_parent'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['neg']['peak_pair_area_heavy'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['neg']['peak_pair_labelp'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['neg']['area_ratio_mask'].insert(10, "Adducts", adduct_flags)
                    self.label_pairs[label]['neg']['peak_pair_corr'].insert(10, "Adducts", adduct_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['neg']['peak_pair_area_parent']["Adducts"] = adduct_flags
                    self.label_pairs[label]['neg']['peak_pair_area_heavy']["Adducts"] = adduct_flags
                    self.label_pairs[label]['neg']['peak_pair_labelp']["Adducts"] = adduct_flags
                    self.label_pairs[label]['neg']['area_ratio_mask']["Adducts"] = adduct_flags
                    self.label_pairs[label]['neg']['peak_pair_corr']["Adducts"] = adduct_flags_corr_sorted
            if polarity not in ['pos', 'neg', 'both']:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    def __flag_adducts_per_label(self, pair_df, area_colnames, polarity, adducts_fnam):
        '''
        Search for, and make a column, for adducts found for each peak pair.
        This function takes an adduct table and will flag adducts
        found for the peak pairs of the polarity specified.
        The adducts are found based on a mass tolerance criterium
        (calculated from ppm_tol), a retention time tolerance criterium
        and lastly a criterium that requires the sum of the adduct peak areas 
        to be smaller than the sum of the parent peak areas.
        The criteria are applied on the parent peaks of each peak pair.
        '''
        # Read adduct table:
        adduct_df = pd.read_csv(adducts_fnam, sep='\t', comment='#')
        # Store adducts in this dict:
        adduct_flag = {i: [] for i in range(len(pair_df))}
        
        # Loop through each pair:
        MW = pair_df['MW_parent'].values
        RT = pair_df['RT_parent'].values
        area_sum = np.sum(pair_df.loc[:, area_colnames].values.astype(float), axis=1)
        for i in range(len(pair_df)):
            # Make a dict with "mass: adduct_name":
            adducts = self.__adduct_expansion(MW[i], adduct_df, polarity)
            # Retention time criterium:
            RT_diff_mask = np.abs(RT[i] - RT) <= self.params['adduct_RT_tol']
            # Smaller area criterium:
            area_mask = area_sum[i] > area_sum

            # Test mass shift criterium, for each possible adduct:
            for MW_adduct in adducts.keys():
                # Mass shift criterium:
                MW_tol = MW_adduct * 1e-6 * self.params['adduct_ppm_tol']
                # Skip adduct if it is M+H or M-H:
                ### Not necessary with area criterium.
                if np.abs(MW_adduct - MW[i]) >= MW_tol:
                    MW_diff_mask = np.abs(MW_adduct - MW) <= MW_tol
                    # Make a mask and store peak pairs:
                    mask = MW_diff_mask & RT_diff_mask & area_mask
                    for idx in np.where(mask)[0]:
                        adduct_flag[idx].append((i, adducts[MW_adduct]))

        # Return list of adduct flags:
        adduct_list = [adduct_flag[idx] if len(adduct_flag[idx])>0 else None for idx in pair_df.index]        
        return(adduct_list)

    def __adduct_expansion(self, MW, adduct_df, polarity):
        '''
        Expand a mass into its adduct masses.
        '''
        MW_adducts = dict()
        for name, mass, charge in zip(adduct_df['Adduct_name'].values, adduct_df['Adduct_mass'].values, adduct_df['Charge'].values):
            if polarity == 'pos' and charge > 0:
                MW_adducts[self.__adduct2mass(MW, name, mass, charge, polarity)] = name
            elif polarity == 'neg' and charge < 0:
                MW_adducts[self.__adduct2mass(MW, name, mass, charge, polarity)] = name

        return(MW_adducts)

    def __adduct2mass(self, MW, adduct_name, adduct_mass, adduct_charge, polarity):
        '''
        For one mass and one adduct calculate the
        M+H or M-H molecular weight.
        '''
        if adduct_name[0].isdigit():
            M_mult = int(adduct_name[0])
            assert(adduct_name[1] == 'M')
        else:
            M_mult = 1
            assert(adduct_name[0] == 'M')

        if polarity == 'pos':
            MW_adduct = M_mult*MW / adduct_charge + adduct_mass - (self.hydrogen_mass - self.electron_mass)
        else:
            MW_adduct = M_mult*MW / -adduct_charge + adduct_mass + (self.hydrogen_mass - self.electron_mass)

        return(MW_adduct)

    # Find isotopes for peak pairs according to polarity:
    def flag_isotopes(self, isotope_set, polarity='both'):   
        for label in sorted(self.labels):        
            if polarity in ['pos', 'both']:
                # Find isotopes:
                isotope_flags = self.__flag_isotopes_per_label(self.label_pairs[label]['pos']['peak_pair_area_parent'], self.area_colnames_pos, isotope_set)
                corr_sort_order_rv = self.label_pairs[label]['pos']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                isotope_flags_corr_sorted = [t[0] for t in sorted(zip(isotope_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['pos']['peak_pair_area_parent'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['pos']['peak_pair_area_heavy'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['pos']['peak_pair_labelp'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['pos']['area_ratio_mask'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['pos']['peak_pair_corr'].insert(10, "Isotopes", isotope_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['pos']['peak_pair_area_parent']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['pos']['peak_pair_area_heavy']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['pos']['peak_pair_labelp']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['pos']['area_ratio_mask']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['pos']['peak_pair_corr']["Isotopes"] = isotope_flags_corr_sorted

            if polarity in ['neg', 'both']:
                # Find isotopes:
                isotope_flags = self.__flag_isotopes_per_label(self.label_pairs[label]['neg']['peak_pair_area_parent'], self.area_colnames_neg, isotope_set)
                corr_sort_order_rv = self.label_pairs[label]['neg']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                isotope_flags_corr_sorted = [t[0] for t in sorted(zip(isotope_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['neg']['peak_pair_area_parent'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['neg']['peak_pair_area_heavy'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['neg']['peak_pair_labelp'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['neg']['area_ratio_mask'].insert(10, "Isotopes", isotope_flags)
                    self.label_pairs[label]['neg']['peak_pair_corr'].insert(10, "Isotopes", isotope_flags)
                except ValueError: # Column already exists
                    self.label_pairs[label]['neg']['peak_pair_area_parent']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['neg']['peak_pair_area_heavy']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['neg']['peak_pair_labelp']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['neg']['area_ratio_mask']["Isotopes"] = isotope_flags
                    self.label_pairs[label]['neg']['peak_pair_corr']["Isotopes"] = isotope_flags
            if polarity not in ['pos', 'neg', 'both']:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    def __flag_isotopes_per_label(self, pair_df, area_colnames, isotope_set):
        '''
        Find and flag isotopes defined in isotope_set.
        The isotopes are found based on a mass tolerance criterium
        (calculated from ppm_tol), a retention time tolerance criterium
        a criterium that requires the sum of the isotope peak areas 
        to be smaller than the sum of the parent peak areas and
        a criterium for minimum correlation between parent and isotope peaks.
        The criteria are applied on the parent peaks of each peak pair.
        '''
        # Store isotopes in this dict:
        isotope_flag = {i: [] for i in range(len(pair_df))}

        # Loop through each pair:
        MW = pair_df['MW_parent'].values
        RT = pair_df['RT_parent'].values
        area = pair_df.loc[:, area_colnames].values.astype(float)
        area_sum = np.sum(area, axis=1)
        for i in range(len(pair_df)):
            # Make a dict with "mass: isotope_name":
            isotopes = {isotope_set[iso_tup]['mass_shift']: ' '.join(iso_tup) for iso_tup in isotope_set}
            # Retention time criterium:
            RT_diff_mask = np.abs(RT[i] - RT) <= self.params['isotope_RT_tol']
            # Smaller area criterium:
            area_mask = area_sum[i] > area_sum
            # Peak area correlation criterium:
            corr_mask = vcorrcoef(area, area[i]) > self.params['isotope_corr_tol']

            # Test mass shift criterium, for each possible isotope:
            for MW_shift_isotope in isotopes.keys():
                MW_isotope = MW_shift_isotope + MW[i]
                # Mass shift criterium:
                MW_tol = MW_isotope * 1e-6 * self.params['isotope_ppm_tol']
                MW_diff_mask = np.abs(MW_isotope - MW) <= MW_tol
                # Make a mask and store peak pairs:
                mask = MW_diff_mask & RT_diff_mask & area_mask & corr_mask
                for idx in np.where(mask)[0]:
                    isotope_flag[idx].append((i, isotopes[MW_shift_isotope]))

        # Return list of isotope flags:
        isotope_list = [isotope_flag[idx] if len(isotope_flag[idx])>0 else None for idx in pair_df.index]

        return(isotope_list)

    # Find blacklisted compounds for peak pairs according to polarity:
    def flag_blacklist(self, blacklist_fnam, polarity='both'):
        # Turn the tsv file into a dictionary:
        blacklist_dict = make_blacklist_dict(blacklist_fnam)
        for label in sorted(self.labels):
            if polarity in ['pos', 'both']:
                # Find blacklisted compounds:
                blacklist_flags = self.__flag_blacklist_per_label(self.label_pairs[label]['pos']['peak_pair_area_parent'], blacklist_dict['pos'])
                corr_sort_order_rv = self.label_pairs[label]['pos']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                blacklist_flags_corr_sorted = [t[0] for t in sorted(zip(blacklist_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['pos']['peak_pair_area_parent'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['pos']['peak_pair_area_heavy'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['pos']['peak_pair_labelp'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['pos']['area_ratio_mask'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['pos']['peak_pair_corr'].insert(10, "Blacklist", blacklist_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['pos']['peak_pair_area_parent']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['pos']['peak_pair_area_heavy']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['pos']['peak_pair_labelp']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['pos']['area_ratio_mask']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['pos']['peak_pair_corr']["Blacklist"] = blacklist_flags_corr_sorted

            if polarity in ['neg', 'both']:
                # Find blacklisted compounds:
                blacklist_flags = self.__flag_blacklist_per_label(self.label_pairs[label]['neg']['peak_pair_area_parent'], blacklist_dict['neg'])
                corr_sort_order_rv = self.label_pairs[label]['neg']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                blacklist_flags_corr_sorted = [t[0] for t in sorted(zip(blacklist_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['neg']['peak_pair_area_parent'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['neg']['peak_pair_area_heavy'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['neg']['peak_pair_labelp'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['neg']['area_ratio_mask'].insert(10, "Blacklist", blacklist_flags)
                    self.label_pairs[label]['neg']['peak_pair_corr'].insert(10, "Blacklist", blacklist_flags)
                except ValueError: # Column already exists
                    self.label_pairs[label]['neg']['peak_pair_area_parent']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['neg']['peak_pair_area_heavy']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['neg']['peak_pair_labelp']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['neg']['area_ratio_mask']["Blacklist"] = blacklist_flags
                    self.label_pairs[label]['neg']['peak_pair_corr']["Blacklist"] = blacklist_flags
            if polarity not in ['pos', 'neg', 'both']:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    def __flag_blacklist_per_label(self, pair_df, blacklist):
        '''
        Find and flag compounds defined in blacklist_dict.
        The compounds are found based on a mass tolerance criterium
        (calculated from ppm_tol) and a retention time tolerance criterium.
        The criteria are applied on the parent peaks of each peak pair.
        '''
        # Store blacklist flags in this dict:
        blacklist_flag = {i: [] for i in range(len(pair_df))}

        MW = pair_df['MW_parent'].values
        RT = pair_df['RT_parent'].values
        blacklist_mask = MW > 0 # dummy mask, all True
        for peak in blacklist:
            MW_i = blacklist[peak]['MW']
            RT_i = blacklist[peak]['RT']
            # Retention time criterium:
            RT_diff_mask = np.abs(RT_i - RT) <= blacklist[peak]['RT_tol']
            # Mass shift criterium:
            MW_tol = MW_i * 1e-6 * blacklist[peak]['MW_ppm_tol']
            MW_diff_mask = np.abs(MW_i - MW) <= MW_tol
            mask = (RT_diff_mask & MW_diff_mask)

            # If found, add description:
            if mask.sum() > 0:
                for idx in np.where(mask)[0]:
                    blacklist_flag[idx].append(blacklist[peak]['Description'])

        # Return list of blacklist flags:
        flag_list = [blacklist_flag[idx] if len(blacklist_flag[idx])>0 else None for idx in pair_df.index]
        
        return(flag_list)

    def flag_labels(self, polarity='both'):
        '''
        For each peak pair for all the labels, use the parent peak MW and RT
        as an identifier and assign to these parent peaks the set of labels
        that they have been observed in.
        '''
        if polarity in ['pos', 'both']:
            # Assign a set of labels to all the parent peaks:
            labelled_peaks = dict()
            for label in sorted(self.labels):
                for peak_pair_id in self.label_pairs[label]['pos']['peak_pair_area_parent']['pair_id'].values:
                    if peak_pair_id[0] in labelled_peaks:
                        labelled_peaks[peak_pair_id[0]].add(label)
                    else:
                        labelled_peaks[peak_pair_id[0]] = {label}

            # Add flags to dataframes:
            for label in sorted(self.labels):
                # Make a list of the labels for each peak pair:
                label_flags = [set2csv(labelled_peaks[peak_pair_id[0]]) for peak_pair_id in self.label_pairs[label]['pos']['peak_pair_area_parent']['pair_id'].values]
                corr_sort_order_rv = self.label_pairs[label]['pos']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                label_flags_corr_sorted = [t[0] for t in sorted(zip(label_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['pos']['peak_pair_area_parent'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['pos']['peak_pair_area_heavy'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['pos']['peak_pair_labelp'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['pos']['area_ratio_mask'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['pos']['peak_pair_corr'].insert(10, "Label_set", label_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['pos']['peak_pair_area_parent']["Label_set"] = label_flags
                    self.label_pairs[label]['pos']['peak_pair_area_heavy']["Label_set"] = label_flags
                    self.label_pairs[label]['pos']['peak_pair_labelp']["Label_set"] = label_flags
                    self.label_pairs[label]['pos']['area_ratio_mask']["Label_set"] = label_flags
                    self.label_pairs[label]['pos']['peak_pair_corr']["Label_set"] = label_flags_corr_sorted

        if polarity in ['neg', 'both']:
            # Assign a set of labels to all the parent peaks:
            labelled_peaks = dict()
            for label in sorted(self.labels):
                for peak_pair_id in self.label_pairs[label]['neg']['peak_pair_area_parent']['pair_id'].values:
                    if peak_pair_id[0] in labelled_peaks:
                        labelled_peaks[peak_pair_id[0]].add(label)
                    else:
                        labelled_peaks[peak_pair_id[0]] = {label}

            # Add flags to dataframes:
            for label in sorted(self.labels):
                # Make a list of the labels for each peak pair:
                label_flags = [set2csv(labelled_peaks[peak_pair_id[0]]) for peak_pair_id in self.label_pairs[label]['neg']['peak_pair_area_parent']['pair_id'].values]
                corr_sort_order_rv = self.label_pairs[label]['neg']['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).index.values
                corr_sort_order = [t[0] for t in sorted(zip(range(len(corr_sort_order_rv)), corr_sort_order_rv), key=lambda x: x[1])]
                label_flags_corr_sorted = [t[0] for t in sorted(zip(label_flags, corr_sort_order), key=lambda x: x[1])]
                try:
                    # Insert them as flags in pair dataframe:
                    self.label_pairs[label]['neg']['peak_pair_area_parent'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['neg']['peak_pair_area_heavy'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['neg']['peak_pair_labelp'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['neg']['area_ratio_mask'].insert(10, "Label_set", label_flags)
                    self.label_pairs[label]['neg']['peak_pair_corr'].insert(10, "Label_set", label_flags_corr_sorted)
                except ValueError: # Column already exists
                    self.label_pairs[label]['neg']['peak_pair_area_parent']["Label_set"] = label_flags
                    self.label_pairs[label]['neg']['peak_pair_area_heavy']["Label_set"] = label_flags
                    self.label_pairs[label]['neg']['peak_pair_labelp']["Label_set"] = label_flags
                    self.label_pairs[label]['neg']['area_ratio_mask']["Label_set"] = label_flags
                    self.label_pairs[label]['neg']['peak_pair_corr']["Label_set"] = label_flags_corr_sorted

        if polarity not in ['pos', 'neg', 'both']:
            raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    # Write the peak pair tables as excel file:
    def write_pairs(self, filename, polarity):
        writer = pd.ExcelWriter('{}.xlsx'.format(filename))
        # Iterate over labels:
        for label in sorted(self.label_pairs.keys()):
            # Assign variables according to polarity:
            if polarity == 'pos':
                peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr = self.label_pairs[label]['pos']['peak_pair_area_parent'], self.label_pairs[label]['pos']['peak_pair_area_heavy'], self.label_pairs[label]['pos']['peak_pair_labelp'], self.label_pairs[label]['pos']['area_ratio_mask'], self.label_pairs[label]['pos']['peak_pair_corr']
            elif polarity == 'neg':
                peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr = self.label_pairs[label]['neg']['peak_pair_area_parent'], self.label_pairs[label]['neg']['peak_pair_area_heavy'], self.label_pairs[label]['neg']['peak_pair_labelp'], self.label_pairs[label]['neg']['area_ratio_mask'], self.label_pairs[label]['neg']['peak_pair_corr']
            else:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))
            # Write to file:
            peak_pair_area_parent.to_excel(writer, sheet_name='parent_area_{}'.format(label))
            peak_pair_corr.to_excel(writer, sheet_name='area_corr_{}'.format(label))
            peak_pair_area_heavy.to_excel(writer, sheet_name='heavy_area_{}'.format(label))
            peak_pair_labelp.to_excel(writer, sheet_name='perc_label_{}'.format(label))
            area_ratio_mask.to_excel(writer, sheet_name='area_ratio_{}'.format(label))

        writer.close()


def set2csv(s):
    return(', '.join(str(si) for si in sorted(s)))

def make_lowercase_dict(obj):
    '''
    Recursive lowercase key/values in dictionary. See:
    https://stackoverflow.com/questions/823030/elegant-pythonic-solution-for-forcing-all-keys-and-values-to-lower-case-in-nest
    '''
    try:
        basestring
    except NameError:
        basestring = (str, bytes)

    if hasattr(obj, 'items'):
        # dictionary
        ret = dict()
        for k, v in obj.items():
            ret[make_lowercase_dict(k)] = make_lowercase_dict(v)
        return(ret)
    elif isinstance(obj, basestring):
        # string
        return obj.lower()
    elif hasattr(obj, '__iter__'):
        # list (or the like)
        ret = list()
        for item in obj:
            ret.append(make_lowercase_dict(item))
        return(ret)
    else:
        # anything else
        return(obj)

def vcorrcoef(X, y):
    '''
    Numpy vectorized correlation coefficient.
    See: https://waterprogramming.wordpress.com/2014/06/13/numpy-vectorized-correlation-coefficient/
    '''
    Xm = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    ym = np.mean(y)
    r_num = np.sum((X-Xm) * (y-ym), axis=1)
    r_den = np.sqrt(np.sum((X-Xm)**2, axis=1) * np.sum((y-ym)**2))
    r = r_num/r_den
    return(r)


def write_filterset(excel_peak_pairs, filename, pair_ppm_tol=2, RT_tol=0.1):
    '''
    Write a filterset for Compound Discoverer to filter
    the peak pairs only.
    '''
    a = '''###
###   This file contains the following filters:
###   
###   Row Filter for Compounds:
###   ------------------------------------
###   OR
###   |  '''

    b = '''###   +--AND
###   |  |  
###   |  +--Molecular Weight is between MW_high and MW_low
###   |  |  
###   |  +--RT [min] is between RT_low and RT_high
###   |  '''

    c = '''###   +--AND
###      |  
###      +--Molecular Weight is between MW_high and MW_low
###      |  
###      +--RT [min] is between RT_low and RT_high
###   ------------------------------------
###   '''

    d = r"""'magellan filter set' 1 'FILENAME'  FiltersetProperties 1  'LastFileName' 'C:\Users\Public\Documents\Thermo\Compound Discoverer 3.0\Common Templates\FilterSets\FILENAME' Filter 'ConsolidatedUnknownCompoundItem' FilterProperties 1  'Filter/DisplayPropertyHint' 'Compounds' 1 NARY_OR N_compounds"""

    e = """NARY_AND 2 isbetween FilterConditionProperties 1  'NamedComparableFilterCondition/DisplayPropertyHint' 'Molecular Weight' property 'System.Double, mscorlib' 'MolecularWeight' constant 'System.Double, mscorlib' 'MW_low' constant 'System.Double, mscorlib' 'MW_high' isbetween FilterConditionProperties 1  'NamedComparableFilterCondition/DisplayPropertyHint' 'RT [min]' property 'System.Double, mscorlib' 'RetentionTime' constant 'System.Double, mscorlib' 'RT_low' constant 'System.Double, mscorlib' 'RT_high'"""


    df = pd.read_excel(excel_peak_pairs)
    mw_list = list(pd.concat([df.loc[:, 'MW_parent'], df.loc[:, 'MW_heavy']]))
    rt_list = list(pd.concat([df.loc[:, 'RT_parent'], df.loc[:, 'RT_heavy']]))
    N_compounds = len(mw_list)

    header = [a]

    tail_1 = d.replace('FILENAME', filename)
    tail_1 = tail_1.replace('N_compounds', str(N_compounds))
    tail = [tail_1]

    for mw, rt in zip(mw_list, rt_list):
        MW_error = mw*pair_ppm_tol*1e-6
        MW_high = round(mw + MW_error, 4)
        MW_low = round(mw - MW_error, 4)

        RT_high = round(rt + RT_tol/2, 2)
        RT_low = round(rt - RT_tol/2, 2)

        s = b.replace('MW_high', str(MW_high))
        s = s.replace('MW_low', str(MW_low))
        s = s.replace('RT_high', str(RT_high))
        s = s.replace('RT_low', str(RT_low))
        header.append(s)

        s = e.replace('MW_high', str(MW_high))
        s = s.replace('MW_low', str(MW_low))
        s = s.replace('RT_high', str(RT_high))
        s = s.replace('RT_low', str(RT_low))
        tail.append(s)

    # Last entry is special:
    s = c.replace('MW_high', str(MW_high))
    s = s.replace('MW_low', str(MW_low))
    s = s.replace('RT_high', str(RT_high))
    s = s.replace('RT_low', str(RT_low))    
    header[-1] = s
    header_str = '\n'.join(header)
    tail_str = ' '.join(tail)

    full_entry = header_str + '\n\n' + tail_str
    with open(filename, 'w') as fh:
        print(full_entry, file=fh)


class Isotopes:
    def __init__(self, IUPAC_atomic_masses, IUPAC_atomic_abundances):
        # Initiate by parsing IUPAC information:
        self.IUPAC_atomic_masses = IUPAC_atomic_masses
        self.IUPAC_atomic_abundances = IUPAC_atomic_abundances
        self.iupac_table = None
        self.isotope_info = None
        self.update_isotope_data() # Fills the iupac_table and isotope_info variable above
        # Make dictionary of isotope abundances 
        self.isotope_abundances = self.__abundance_dict()
        self.relevant_isotopes = ['[2]H', '[13]C', '[15]N', '[17]O', '[18]O', '[33]S', '[34]S', '[36]S']
        self.iso_set = None
        
    def find_iso_set(self, iso_names=None, min_abs=1e-7):
        '''Function wrapper, to avoid global variables.'''
        if iso_names is None:
            iso_names = self.relevant_isotopes
        iso_set = dict()
        def new_iso(iso_set, combo=()):
            '''
            Recursive function to find all isotope combinations
            below a certain minimum "expected abundance", based on
            input abundance.
            The "expected abundance" should be understood simply as
            the product of the isotope abundances for a given combination
            of isotopes e.g. ('13C', '13C', '15N') has an expected abundance
            of ~4.477e-07 based on natural isotope abundance.
            Be aware that this is __not__ the same as saying that 13C2 15N of
            a molecule has abundance 4.477e-07 of its parent. This depends
            on the number of carbon and nitrogen in the molecule which
            is unknown a priori.
            Is the formula known, the true expected abundance can be
            calculated using the binomial pmf:
            https://en.wikipedia.org/wiki/Binomial_distribution#Probability_mass_function
            '''
            for iso_name in iso_names:
                combo += (iso_name, )

                # Already in set:
                if combo in iso_set:
                    combo = combo[0:-1]

                # First element:
                elif len(combo) == 1 and self.isotope_abundances[iso_name] >= min_abs:
                    iso_set[combo] = self.isotope_abundances[iso_name]
                    new_iso(iso_set, combo=combo) # Recursion, add more 

                # Move to next isotope:
                elif len(combo) == 1 and self.isotope_abundances[iso_name] < min_abs:
                    combo = combo[0:-1]

                # Additional elements:
                elif combo not in iso_set and (iso_set[combo[0:-1]] * self.isotope_abundances[iso_name]) >= min_abs:
                    iso_set[combo] = iso_set[combo[0:-1]] * self.isotope_abundances[iso_name]
                    new_iso(iso_set, combo=combo) # Recursion, add more 

                # Move to next isotope:
                elif combo not in iso_set and (iso_set[combo[0:-1]] * self.isotope_abundances[iso_name]) < min_abs:
                    combo = combo[0:-1]

                # All outcomes have been expanded with "elif"
                # to increase clarity, fallback on "else" should not happen:
                else:
                    assert(False) # should not happen

            if len(combo) > 0:
                new_iso(iso_set, combo=combo[0:-1]) # Peel off one layer, then add more
            return(iso_set)

        # Run recursion:
        iso_set = new_iso(iso_set, combo=())
        # Remove duplicates:
        iso_set_sorted = dict()
        for combo in iso_set.keys():
            sorted_combo = tuple(sorted(combo))
            iso_set_sorted[sorted_combo] = iso_set[combo]
        
        # Add mass shift to dictionary:
        for combo, abundance in iso_set_sorted.items():
            iso_set_sorted[combo] = {'abundance': abundance, 'mass_shift': self.isotopes2mass_shift(' '.join(combo))}
        
        self.iso_set = iso_set_sorted
        return(iso_set_sorted)

    def __abundance_dict(self):
        abundances = dict()
        for element in self.isotope_info.keys():
            for nominal_mass in self.isotope_info[element].keys():
                if nominal_mass != 'most_abundant':
                    nuclide = '[{}]{}'.format(nominal_mass, element)
                    abundances[nuclide] = self.isotope_info[element][nominal_mass]['abundance']
        return(abundances)

    def isotopes2mass_shift(self, isotope_str):
        '''
        Convert an isotope string to a mass shift.
        E.g. "[13]C3 [15]N" will return the mass shift from 3x 13C carbon and 1x 15N nitrogen.
        '''
        labelled_mass = self.formula2mass(isotope_str)
        isotope_str_unlabelled = list()
        for element in isotope_str.split():
            if '[' == element[0]:
                end = element.index(']')
                element = element[(end+1):]
            isotope_str_unlabelled.append(element)
        unlabelled_mass = self.formula2mass(' '.join(isotope_str_unlabelled))

        return(labelled_mass - unlabelled_mass)

    def formula2mass(self, formula):
        '''
        Given a whitespace separated chemical formula, calculate the exact mass.
        Isotopes deviating from the naturally most abundant are indicated by 
        their nominal (integer) mass enclodes by brackets.
        E.g. 13C labelled glucose is: [13]C6 H12 O6 
        A mix of labelled and unlabelled is written separated e.g.
        glucose with half of carbons labelled is: C3 [13]C3 H12 O6 
        '''
        total_mass = 0
        for element in formula.split():
            symbol_not_found = True
            for symbol in self.isotope_info.keys():
                # No bracket means this is the most abundant isotope:
                if symbol in element and '[' != element[0]:
                    symbol_not_found = False
                    split_element = element.split(symbol)
                    if split_element[0] != '':
                        raise Exception('Wrong formatting of formula string: {}.\nEither start with bracket and the nominal mass, or an element abbreviation.'.format(element))
                        
                    if split_element[1] != '':
                        try:
                            N = float(split_element[1])
                        except ValueError: # symbol not done
                            continue
                    else:
                        N = 1
                    total_mass += N * self.isotope_info[symbol]['most_abundant']['mass']

                # Bracket means we need to parse the isotope nominal mass:
                elif symbol in element and '[' == element[0]:
                    symbol_not_found = False
                    try:
                        end = element.index(']')
                        nominal_mass = int(element[1:end])
                        split_element = element[(end+1):].split(symbol)
                    except:
                         raise Exception('Could not read the input formula: {}\nPlease check formating.'.format(element))
                    if split_element[1] != '':
                        try:
                            N = float(split_element[1])
                        except ValueError: # symbol not done
                            continue
                    else:
                        N = 1
                    try:
                        total_mass += N * self.isotope_info[symbol][nominal_mass]['mass']
                    except KeyError:
                        raise Exception('Could not find the nominal mass of {} for this element ({}), please check the input formula: {}'.format(nominal_mass, symbol, element))

            if symbol_not_found:
                raise Exception('Did not find an element symbol in this part of the input formula: {}\nUse standard element abbreviations (first letter uppercase, second letter lower) and separate each element by a space.'.format(element))
                
        return(total_mass)

    def update_isotope_data(self):
        self.iupac_table, self.isotope_info = self.__read_iupac(self.IUPAC_atomic_masses, self.IUPAC_atomic_abundances)

    def __read_iupac(self, IUPAC_atomic_masses, IUPAC_atomic_abundances):
        '''
        Parse csv and html tables from IUPAC:
        https://iupac.org/isotopes-matter/
        Containing information about atomic masses and
        isotopic abundance.
        Then return this as a dataframe and a dictionary.
        '''
        # The csv file with atomic masses contains
        # html links with the year when the measurement was made
        # therefore make a parser to extract this year:
        from html.parser import HTMLParser
        class MyHTMLParser(HTMLParser):
            def __init__(self):
                super().__init__()
                self.data_extract = list()
            def handle_data(self, data):
                self.data_extract.append(int(data))
        html_parser = MyHTMLParser()
        
        # Read the IUPAC masses into a dataframe:
        iupac_masses = pd.read_csv(IUPAC_atomic_masses, sep=',', comment='#')
        # Extract the year of measurement:
        for html_link in iupac_masses['Year/link'].values:
            html_parser.feed(html_link)
        iupac_masses['year'] = html_parser.data_extract
        iupac_masses = iupac_masses.drop(['Year/link', 'uncertainty'], axis=1)
        # Extract data row from the most recent year of measurement,
        # then sort by mass (smallest to largest):
        iupac_masses_agg = iupac_masses.iloc[iupac_masses.groupby('nuclide').agg({'year': 'idxmax'}).year].sort_values(by=['mass']).reset_index(drop=True)
        iupac_masses_agg = iupac_masses_agg.drop(['year'], axis=1)

        # Read the IUPAC abundances into a dataframe:
        iupac_abundances = pd.read_html(IUPAC_atomic_abundances, encoding='utf-8')[0]
        Z_number = iupac_abundances['Z'].str.isnumeric().fillna(False)
        iupac_abundances = iupac_abundances[Z_number]
        iupac_abundances['nuclide'] = [A+E for A, E in zip(iupac_abundances['A'].values, iupac_abundances['E'].values)]

        # Abundance data comes as either:
        # 1) a range, 2) a point estimate, or 3) no info
        # Calculate a point estimate for all, dump the ones with no info:
        mean_abundances = list()
        for ric in iupac_abundances['Representative isotopic composition'].values:
            # Remove whitespace:
            ric = ''.join(ric.split())
            try:
                # Ignore warnings:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="'float' object is not callable; perhaps you missed a comma?")
                    ric_range = eval(ric)
                    ric_point_est = sum(ric_range)/2
            except TypeError:
                ric_point_est = float(ric.split('(')[0])
            except SyntaxError:
                ric_point_est = np.nan
            mean_abundances.append(ric_point_est)
        iupac_abundances['abundance'] = mean_abundances
        abs_mask = ~iupac_abundances['abundance'].isna()
        iupac_abundances = iupac_abundances[abs_mask]
        iupac_abundances = iupac_abundances.drop(['Notes', 'Representative isotopic composition'], axis=1)
        
        # Merge masses with abundances:
        iupac_isotopes = iupac_abundances.merge(iupac_masses_agg, on='nuclide')
        iupac_isotopes['A'] = [int(val) for val in iupac_isotopes['A'].values]
        iupac_isotopes['Z'] = [int(val) for val in iupac_isotopes['Z'].values]

        # Convert this into a dictionary usefull for
        # converting a chemical formula into a mass:
        isotope_info = dict()
        for element, nominal_mass, mass, abundance in zip(iupac_isotopes['E'].values, iupac_isotopes['A'].values, iupac_isotopes['mass'].values, iupac_isotopes['abundance'].values):
            if element not in isotope_info:
                isotope_info[element] = {nominal_mass: {'mass': mass, 'abundance': abundance}}
            else:
                isotope_info[element][nominal_mass] = {'mass': mass, 'abundance': abundance}

        # Insert a special key for the most abundant isotope:
        for element in isotope_info.keys():
            abundance_sorter = [(isotope_info[element][nominal_mass]['abundance'], nominal_mass) for nominal_mass in isotope_info[element].keys()]
            most_abundant = sorted(abundance_sorter, key=lambda x: x[0], reverse=True)
            isotope_info[element]['most_abundant'] = isotope_info[element][most_abundant[0][1]]  

        return(iupac_isotopes, isotope_info)

def make_blacklist_dict(blacklist_fnam):
    blacklist_df = pd.read_csv(blacklist_fnam, sep='\t')
    blacklist_dict = {'pos': {}, 'neg': {}}
    blacklist_cols = blacklist_df.columns
    for i in range(len(blacklist_df)):
        peak_id = (blacklist_df.loc[i, 'MW'], blacklist_df.loc[i, 'RT'])
        assert(peak_id not in blacklist_dict)
        if blacklist_df.loc[i, 'Polarity'] == 'pos':
            blacklist_dict['pos'][peak_id] = {col: blacklist_df.loc[i, col] for col in blacklist_cols}
        elif blacklist_df.loc[i, 'Polarity'] == 'neg':
            blacklist_dict['neg'][peak_id] = {col: blacklist_df.loc[i, col] for col in blacklist_cols}
        else:
            raise Exception('Polarity "{}" not recognized'.format(blacklist_df.loc[i, 'Polarity']))
    return(blacklist_dict)


