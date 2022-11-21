import pandas as pd
import numpy as np
import copy

class PeakData:
    import pandas as _pd
    import numpy as _np
    import random as _rd
#    from IPython.display import display, HTML
    import copy as _cp

    def __init__(self, name, sample_names, sample_label, min_label, params):
        self.name = name
        self.sample_names = sample_names # Dictionary mapping the area column (sample) names to a more descriptive sample name
        self.sample_label = sample_label # Which sample names are labelled with what label e.g. ?
        self.min_label = min_label       # Minimum number of labeled samples where the peak pair passes area ratio criterium
        self.params = params
        # Dummy for peak data:
        self.peak_data_pos = None
        self.peak_data_neg = None
        # Dummy for column names containing peak area:
        self.area_colnames_pos = None
        self.area_colnames_neg = None

        # For each label make a dictionary entry with all the
        # tables related to the peak pair from that label:
        self.labels = list(sample_label.keys())
        self.label_peaks = dict()
        # Iterate over labels:
        for label in sorted(self.labels):
            self.label_peaks[label] = dict()
            self.label_peaks[label]['label_colnames'] = sample_label[label]
            for polarity in ['pos', 'neg']:
                self.label_peaks[label][polarity] = dict()
                self.label_peaks[label][polarity]['peak_pair_area_parent'] = None # Parent peak area
                self.label_peaks[label][polarity]['peak_pair_area_heavy'] = None # Labelled peak area
                self.label_peaks[label][polarity]['peak_pair_labelp'] = None # Labelling percent
                self.label_peaks[label][polarity]['area_ratio_mask'] = None # Filter based on self.params e.g. mass different, retention time difference etc.
                self.label_peaks[label][polarity]['peak_pair_corr'] = None # Matrix with accross sample correlation coefficients between peak areas

    # Read peak data according to polarity:
    def read_peaks(self, datafile, polarity):
        if polarity == 'pos':
            self.peak_data_pos = self.__peak_reader(datafile, polarity)
        elif polarity == 'neg':
            self.peak_data_neg = self.__peak_reader(datafile, polarity)
        else:
            raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

    # Internal function to read peaks in either polarity:
    def __peak_reader(self, datafile, polarity):
        csv = self._pd.read_excel(datafile)
        csv.columns = [n.lower().replace(' ', '_') for n in csv.columns]
        # Assert that weight and retention time can be used as unique identifier:
        wr_list = [(w, r) for w, r, in zip(csv['molecular_weight'], csv['rt_[min]'])]
        assert(len(wr_list) == len(set(wr_list)))
        # Insert (MW, RT) as peak ID into the dataframe:
        csv.insert(0, 'MW', csv['molecular_weight'])
        csv.insert(0, 'RT', csv['rt_[min]'])
        csv.insert(0, 'peak_id', self._pd.Series(wr_list, index=csv.index))
        # Extract the columns with area information:
        col_sele = ['peak_id', 'RT', 'MW', 'name']
        col_sele.extend([col for col in csv.columns if 'area' in col])
        peak_data = csv.loc[:, col_sele]

        # Get sample names from sample description:
        sample_names = list(self.sample_names.keys())
        colname_list = []
        area_colnames = []
        for col in peak_data.columns:
            if 'area:' in col:
                colname = False
                for name in sample_names:
                    if colname is False and name.lower() in col:
                        colname = self.sample_names[name]
                    elif colname is True and name.lower() in col:
                        raise Exception('Column "{}" maps to multiple names in sample description. Please correct.'.format(col))
                if colname is False:
                        raise Exception('Column name "{}" not found. Check sample description.'.format(col))
                colname_list.append(colname)
                area_colnames.append(colname)
            else:
                colname_list.append(col)
        peak_data.columns = colname_list
        
        if polarity == 'pos':
            self.area_colnames_pos = area_colnames
        elif polarity == 'neg':
            self.area_colnames_neg = area_colnames
        
        # Filter peak dataframe:
        peak_data = self.__run_all_peak_filters(peak_data, area_colnames, self.params)
        peak_data = peak_data.sort_values(by='MW', ascending=True)
        peak_data.reset_index(drop=True, inplace=True)

        return(peak_data)

    # Apply all filters on the peak data:
    def __run_all_peak_filters(self, peak_data, area_colnames, params):
        '''
        Function to run sequential filtering on peaks.
        Currently, only has one filter, but more could be added.
        '''
        count_before = len(peak_data)
        peak_data = self.min_area_peak_filter(peak_data, area_colnames, self.params['min_area'])
        count_after = len(peak_data)
        print('Filtered {} peaks out. {} peaks left.'.format(count_before - count_after, count_after))
        return(peak_data)

    # Filter peaks with area smaller than the cutoff:
    def min_area_peak_filter(self, peak_data, area_colnames, min_area):
        '''
        Enforce a minimum peak area.
        For a given peak the peak area at least one sample has to be above this value.
        '''
        area_max = peak_data.loc[:, area_colnames].max(axis=1)
        mask = (area_max > min_area)
        peak_data_filtered = peak_data[mask]
        return(peak_data_filtered)

    # Find peak pairs according to polarity:
    def find_pairs(self, polarity):
        for label in sorted(self.labels):            
            if polarity == 'pos':
                self.label_peaks[label]['pos']['peak_pair_area_parent'], self.label_peaks[label]['pos']['peak_pair_area_heavy'], self.label_peaks[label]['pos']['peak_pair_labelp'], self.label_peaks[label]['pos']['area_ratio_mask'], self.label_peaks[label]['pos']['peak_pair_corr'] = self.__find_peak_pairs_per_label(self.peak_data_pos, self.area_colnames_pos, polarity, label)
            elif polarity == 'neg':
                self.label_peaks[label]['neg']['peak_pair_area_parent'], self.label_peaks[label]['neg']['peak_pair_area_heavy'], self.label_peaks[label]['neg']['peak_pair_labelp'], self.label_peaks[label]['neg']['area_ratio_mask'], self.label_peaks[label]['neg']['peak_pair_corr'] = self.__find_peak_pairs_per_label(self.peak_data_neg, self.area_colnames_neg, polarity, label)
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
                          'polarity', 'label', 'name']
        pair_columns = self._cp.deepcopy(pair_info_mask)
        # Add area to the peak pair table:
        pair_columns.extend(area_colnames)
        # Make a table for the area of parent/heavy compound:
        peak_pair_area_parent = self._pd.DataFrame(columns = pair_columns)
        peak_pair_area_heavy = self._pd.DataFrame(columns = pair_columns)

        peak_pair_set = set()
        ppm_cutoff = self.params['MW_shift_ppm_tol']
        RT_tol = self.params['RT_tol']
        MW = peak_data.loc[:, 'MW'].values
        RT = peak_data.loc[:, 'RT'].values
        peak_id = peak_data.loc[:, 'peak_id'].values
        MW_shift = self.params['MW_shift'][label] # The mass shift is defined by the label
        # Iterate over peaks to match on MW/RT criteria:
        for i in range(len(peak_data)):
            # Retention time criterium:
            RT_diff_mask = self._np.abs(RT[i] - RT) <= RT_tol
            # Mass shift criterium:
            MW_tol = MW[i] * 1e-6 * ppm_cutoff
            MW_diff_high = MW_shift + MW_tol
            MW_diff_low = MW_shift - MW_tol
            MW_diff_higer_mask = ((MW - MW[i]) <= MW_diff_high) & ((MW - MW[i]) >= MW_diff_low)

            # Make a mask and store peak pairs:
            mask = (RT_diff_mask & MW_diff_higer_mask)
            if mask.sum() > 0:
                for j in self._np.where(mask)[0]:
                    # Because df is sorted by MW,
                    # i is parent, j is heavy
                    CD_name = peak_data.loc[i, 'name']
                    pair_id = (peak_id[i], peak_id[j], polarity, label)
                    assert(pair_id not in peak_pair_set)
                    peak_pair_set.add(pair_id)

                    # Add data to each peak pair table:
                    tmp_df = self._pd.DataFrame(columns = pair_columns)
                    tmp_df.loc[i, pair_info_mask] = self._np.array([pair_id, MW[i], RT[i], MW[j], RT[j], polarity, label, CD_name], dtype=object)
                    tmp_df.loc[i, area_colnames] = peak_data.loc[i, area_colnames]
                    peak_pair_area_parent = peak_pair_area_parent.append(tmp_df, ignore_index=True)
                    tmp_df = self._pd.DataFrame(columns = pair_columns)
                    tmp_df.loc[j, area_colnames] = peak_data.loc[j, area_colnames]
                    tmp_df.loc[j, pair_info_mask] = self._np.array([pair_id, MW[i], RT[i], MW[j], RT[j], polarity, label, CD_name], dtype=object)
                    peak_pair_area_heavy = peak_pair_area_heavy.append(tmp_df, ignore_index=True)

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
                               (peak_pair_labelp <= area_ratio_cutoff[1]))) # ratio must be below max cutoff

        # Add pair info:
        pair_info_df = peak_pair_area_heavy.loc[:, pair_info_mask]
        peak_pair_labelp = self._pd.concat([pair_info_df, peak_pair_labelp], axis=1, sort=False)
        area_ratio_mask = self._pd.concat([pair_info_df, area_ratio_mask], axis=1, sort=False)

        # Find and drop rows with insufficient number of samples passing the area ratio criterium:
        area_ratio_drop = list()
        for i in range(len(peak_pair_area_parent)):
            cols_with_label = area_ratio_mask.columns.isin((self.label_peaks[label]['label_colnames']))
            if sum(area_ratio_mask.loc[i, cols_with_label]) < self.min_label:
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
        peak_pair_corr = self._pd.concat([pair_info_df, peak_pair_corr], axis=1, sort=False)

        return(peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr)

    # Filter all peaks not in the intersection of
    # a given set of labels:
    def intersection_pairs(self, labels, polarity):
        if polarity == 'pos':
            area_colnames = self.area_colnames_pos
        elif polarity == 'neg':
            area_colnames = self.area_colnames_neg
        else:
            raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))

        pair_info_mask = ['pair_id', 'MW_parent', 'RT_parent', 'MW_heavy', 'RT_heavy',
                          'polarity', 'label', 'name']
        # Find the intersection based on parent MW and RT:
        s_intersection = set([pi[0] for pi in self.label_peaks[labels[0]][polarity]['peak_pair_area_parent']['pair_id'].values])
        for label in labels:
            s = set([pi[0] for pi in self.label_peaks[label][polarity]['peak_pair_area_parent']['pair_id'].values])
            s_intersection = s_intersection.intersection(s)
        
        # Update peak pair dataframe for the chosen labels:
        for label in labels:
            row_mask = [pi[0] in s_intersection for pi in self.label_peaks[label][polarity]['peak_pair_area_parent']['pair_id'].values]
            
            self.label_peaks[label][polarity]['peak_pair_area_parent'] = self.label_peaks[label][polarity]['peak_pair_area_parent'].loc[row_mask, :].reset_index(drop=True, inplace=False)
            self.label_peaks[label][polarity]['peak_pair_area_heavy'] = self.label_peaks[label][polarity]['peak_pair_area_heavy'].loc[row_mask, :].reset_index(drop=True, inplace=False)
            self.label_peaks[label][polarity]['peak_pair_labelp'] = self.label_peaks[label][polarity]['peak_pair_labelp'].loc[row_mask, :].reset_index(drop=True, inplace=False)
            self.label_peaks[label][polarity]['area_ratio_mask'] = self.label_peaks[label][polarity]['area_ratio_mask'].loc[row_mask, :].reset_index(drop=True, inplace=False)
            
            # Redo the dataframe with peak pair area correlations:
            a = self.label_peaks[label][polarity]['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).loc[:, area_colnames]
            a = a.astype(float)
            peak_pair_corr = a.T.corr().abs()
            pair_info_df = self.label_peaks[label][polarity]['peak_pair_area_parent'].sort_values(by='RT_parent', ascending=True).loc[:, pair_info_mask]
            self.label_peaks[label][polarity]['peak_pair_corr'] = self._pd.concat([pair_info_df, peak_pair_corr], axis=1, sort=False)

    # Write the peak pair tables as excel file:
    def write_pairs(self, filename, polarity):
        writer = self._pd.ExcelWriter('{}.xlsx'.format(filename))
        # Iterate over labels:
        for label in sorted(self.label_peaks.keys()):
            # Assign variables according to polarity:
            if polarity == 'pos':
                peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr = self.label_peaks[label]['pos']['peak_pair_area_parent'], self.label_peaks[label]['pos']['peak_pair_area_heavy'], self.label_peaks[label]['pos']['peak_pair_labelp'], self.label_peaks[label]['pos']['area_ratio_mask'], self.label_peaks[label]['pos']['peak_pair_corr']
            elif polarity == 'neg':
                peak_pair_area_parent, peak_pair_area_heavy, peak_pair_labelp, area_ratio_mask, peak_pair_corr = self.label_peaks[label]['neg']['peak_pair_area_parent'], self.label_peaks[label]['neg']['peak_pair_area_heavy'], self.label_peaks[label]['neg']['peak_pair_labelp'], self.label_peaks[label]['neg']['area_ratio_mask'], self.label_peaks[label]['neg']['peak_pair_corr']
            else:
                raise Exception('The polarity "{}" could not be recognized, not pos/neg.'.format(polarity))
            # Write to file:
            peak_pair_area_parent.to_excel(writer, sheet_name='parent_area_{}'.format(label))
            peak_pair_corr.to_excel(writer, sheet_name='area_corr_{}'.format(label))
            peak_pair_area_heavy.to_excel(writer, sheet_name='heavy_area_{}'.format(label))
            peak_pair_labelp.to_excel(writer, sheet_name='perc_label_{}'.format(label))
            area_ratio_mask.to_excel(writer, sheet_name='area_ratio_{}'.format(label))

        writer.close()


def pick_ratio(peak_data, area_colnames, labeled_columns, known_fnam, label_str, params):
    label_df = pd.read_csv(known_fnam, sep='\t')
    compounds = sorted(list(set(label_df['Name'].values)))
    data_idx = [c+'_'+i for c in compounds for i in label_str]
    data_idx = [n + "_" + f for n, f in zip(label_df['Name'], label_df['Label'])]
    perc_lab = pd.DataFrame(columns=area_colnames, index=data_idx)
    for compound in compounds:
        dp_df = label_df[label_df['Name'] == compound]
        area = {label_str[0]:[0], label_str[1]:[0]}

        for index, row in dp_df.iterrows():
            mask = df_filter(peak_data, area_colnames, row['RT'], row['Mass'], params)
            if sum(mask) == 1:
                area[row['Label']] = peak_data[mask][area_colnames].values[0]

        if sum(area[label_str[0]]) > 0 and sum(area[label_str[1]]) > 0:
            area_sum = area[label_str[0]] + area[label_str[1]]
            perc_m = area[label_str[0]] / area_sum
            perc_m4 = area[label_str[1]] / area_sum

            idx = compound+'_'+label_str[0]
            perc_lab.loc[idx] = perc_m

            idx = compound+'_'+label_str[1]
            perc_lab.loc[idx] = perc_m4

    #print(perc_lab)
    x_labeled = perc_lab.loc[:, perc_lab.columns.isin((labeled_columns))].dropna().astype('float64')
    desc = x_labeled.T.describe(include='all')

    return(desc)


def df_filter(df, area_colnames, RT, Mass, params):
    Mass_error = params['MW_shift_ppm_tol'] * Mass * 1e-6
    RT_mask = ((RT + params['RT_tol']) > df['RT']) & ((RT - params['RT_tol']) < df['RT'])
    MW_mask = ((Mass + Mass_error) > df['MW']) & ((Mass - Mass_error) < df['MW'])
    mask = RT_mask & MW_mask
    mask_area_max = df[area_colnames].sum(1) == df[mask][area_colnames].sum(1).max()
    mask = mask & mask_area_max

    return(mask)


def write_filterset(filter_data, filename):
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


    df = pd.read_excel(filter_data)
    mw_list = list(pd.concat([df.loc[:, 'MW_parent'], df.loc[:, 'MW_heavy']]))
    rt_list = list(pd.concat([df.loc[:, 'RT_parent'], df.loc[:, 'RT_heavy']]))
    N_compounds = len(mw_list)

    header = [a]

    tail_1 = d.replace('FILENAME', filename)
    tail_1 = tail_1.replace('N_compounds', str(N_compounds))
    tail = [tail_1]

    for mw, rt in zip(mw_list, rt_list):
        MW_error = mw*params['MW_shift_ppm_tol']*1e-6
        MW_high = round(mw + MW_error, 4)
        MW_low = round(mw - MW_error, 4)

        RT_high = round(rt + params['RT_tol']/2, 2)
        RT_low = round(rt - params['RT_tol']/2, 2)

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

