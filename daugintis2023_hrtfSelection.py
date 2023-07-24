from pathlib import Path
import os
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patches as mpatches

my_cmap = ['firebrick', 'mistyrose', 'white', 'lightgreen', 'darkgreen']
handles = [plt.Rectangle((0, 0), 0, 0, facecolor='white',
                         edgecolor='white'),
           plt.Rectangle((0, 0), 0, 0, facecolor=my_cmap[1],
                         edgecolor='black'),
           plt.Rectangle((0, 0), 0, 0, facecolor=my_cmap[1],
                         edgecolor=my_cmap[0], hatch='///'),
           plt.Rectangle((0, 0), 0, 0, facecolor=my_cmap[3],
                         edgecolor='black'),
           plt.Rectangle((0, 0), 0, 0, facecolor=my_cmap[3],
                         edgecolor=my_cmap[4], hatch='\\\\\\'),
           plt.Rectangle((0, 0), 0, 0, facecolor=my_cmap[2],
                         edgecolor='black'), ]
my_labels = ['Worst', 'Bad', 'Individual', 'Good', 'Best']

def select_HRTFs(plot_figures=False):
    """
    Selects the best and the worst HRTF  based on the dir_err_table.csv file in the directory
    :param plot_figures: if True than the error distributions for each subject with all HRTF
     will be plotted and saved in the folder
    :return: Dataframe with the best and the worst HRTF for each subject.
    """
    dir_err_files = Path('model_predictions').glob('dir_err_table_*.csv')
    
    selected_hrtf_list = pd.DataFrame(columns=['id', 'goodHRTF', 'badHRTF'])
    if plot_figures:
        good_err_all = pd.DataFrame(columns=['template_HRTF', 'target_HRTF'])
        bad_err_all = pd.DataFrame(columns=['template_HRTF', 'target_HRTF'])
        if not os.path.exists('selection_figures'):
            os.mkdir('selection_figures')
    
    for dir_err_file in dir_err_files:
        print(f'%%%%%%%%%%% Reading file: {dir_err_file} %%%%%%%%%%%')
        dir_errors = pd.read_csv(dir_err_file, sep=',', engine='c', header=0)
        dir_errors = dir_errors[(np.abs(dir_errors.pol_target) <= 11.5) & (np.abs(dir_errors.lat_target) <= 30)]

        template_hrtfs = dir_errors["template_HRTF"].unique()
        target_hrtfs = dir_errors["target_HRTF"].unique()
        for template in template_hrtfs:
            print(f'%%%% Template: {template} %%%%%')

            good_errors, bad_errors = classify_hrtf_distributions(template, target_hrtfs, dir_errors)
            best_hrtf, worst_hrtf = select_best_worst_hrtf(good_errors, bad_errors)

            selected_hrtf_list.loc[len(selected_hrtf_list)] = [template, best_hrtf, worst_hrtf]

            if plot_figures:
                good_err_all = pd.concat([good_err_all, good_errors])
                bad_err_all = pd.concat([bad_err_all, bad_errors])
                plot_hrtf_distributions(template, target_hrtfs,
                                        good_errors, best_hrtf,
                                        bad_errors, worst_hrtf,
                                        dir_errors,
                                        save_figure=True)


    selected_hrtf_list.to_csv('selected_HRTFs.csv')

    if plot_figures:
        sel_m = selection_matrix(selected_hrtf_list, template_hrtfs, target_hrtfs, good_err_all, bad_err_all)
        # anonim_ticks = np.arange(5,sel_m.columns.size+1,5)
        # sns.set(font_scale=1)
        sns.set_theme(style="whitegrid")
        # for s, t in zip([st_ind, st_avg], ['ind_selection', 'median_selection']):
        cm = 1/2.54
        plt.rcParams.update({'font.size': 10})
        fig = plt.figure(figsize=[7.8*cm, 7.8*cm])
        ax = sns.heatmap(data=sel_m, linewidths=.5, cbar=False, mask=(sel_m == 0), cmap=my_cmap, vmin=-2, vmax=2)
        ax.invert_yaxis()
        # ax.set_xticks(anonim_ticks-.5,anonim_ticks)
        # ax.set_yticks(anonim_ticks-.5,anonim_ticks)
        ax.tick_params(axis='both', labelsize='small', pad=-5, labelrotation=0)
        ax.set_xlabel("Subject")
        ax.set_ylabel("Target HRTF")
        ax.grid(False)
        pch = []
        for n, cl in enumerate(my_cmap):
            if n == 2:
                continue
            pch.append(mpatches.Patch(color=cl, label=my_labels[n]))
        pch[0], pch[1] = pch[1], pch[0]
        fig.legend(handles=pch, ncol=4, loc='upper center', fontsize='small', columnspacing=1, bbox_to_anchor=(0.48, 0.97), frameon=False, handletextpad=.4)
        plt.savefig(Path('selection_figures')/Path('selection_matrix.pdf'), format='pdf', bbox_inches='tight')

    return selected_hrtf_list



def classify_hrtf_distributions(subject: str, target_hrtfs: np.ndarray, errors_lim_dirs: pd.DataFrame):
    """
    Classifies target HRTFs for a subject into good and bad based on Shapiro-Wilk normality test.
    :param subject: subject id
    :param target_hrtfs: array of target HRTF ids
    :param errors_lim_dirs: Dataframe with directional errors, limited to the desired directions
    :return: tuple of dataframes with errors of good target HRTFs and bad target HRTFs for the subject
    """
    good_errors = pd.DataFrame(columns=['template_HRTF', 'target_HRTF', 'rmsP', 'querr_3rd_quartile'])
    bad_errors = pd.DataFrame(columns=['template_HRTF', 'target_HRTF', 'querr_median', 'querr_3rd_quartile', 'rmsP'])
    n_good, n_bad = 0, 0  # Counters for good/bad HRTF tables
    for target in target_hrtfs:
        if target != subject:
            # Filtering the directional errors for the specific subject and target
            target_errors = errors_lim_dirs[(errors_lim_dirs['template_HRTF'] == subject) &
                                            (errors_lim_dirs['target_HRTF'] == target)]
            shapiro_test = stats.shapiro(target_errors.rmsP)
            # Good HRTFs
            if shapiro_test.pvalue > 0.05:
                # Aggregated global rms polar error
                rms_p_global = rms_p(target_errors['rmsP'])
                # 3rd quartile of quadrant error
                q_err_3q = target_errors.querr.quantile(0.75)
                good_errors.loc[n_good] = [subject, target, rms_p_global, q_err_3q]
                n_good += 1
            # Bad HRTFs
            else:
                # Median quadrant error
                q_err_median = target_errors.querr.median()
                # 3rd quartile of qadrant error
                q_err_3q = target_errors.querr.quantile(0.75)
                # Global rms polar error
                rms_p_global = rms_p(target_errors['rmsP'])
                bad_errors.loc[n_bad] = [subject, target, q_err_median, q_err_3q, rms_p_global]
                n_bad += 1
    return good_errors, bad_errors


def select_best_worst_hrtf(good_errors_df: pd.DataFrame, bad_errors_df: pd.DataFrame):
    """
    Function that selects the 'best' and the 'worst' HRTFs from respective
    error tables for the subject.
    :param good_errors_df: Error distributions for target HRTFs
        classed as 'good'. Dataframe must contain columns with 'target_HRTF',
        'rmsP', and 'querr_3rd_quartile' columns.
    :param bad_errors_df: Dataframe of error distributions for target HRTFs
        classed as 'bad'. Dataframe must contain columns with 'target_HRTF',
        'querr_median', 'querr_3rd_quartile', and 'rmsP' columns.
    :return: tuple with good and bad HRTF labels
    """

    # Finding the 'best' HRTF
    n = 0
    good_hrtf = set()
    while len(good_hrtf) < 1:
        n += 1
        if not good_errors_df.empty:
            # Selecting to n top HRTFs based on PE
            n_good_pe = good_errors_df.nsmallest(n, 'rmsP', keep='all')
            # Selecting to n top HRTFs based on QE
            n_good_qe = good_errors_df.nsmallest(n, 'querr_3rd_quartile', keep='all')
            # Looking for common HRTFs between top n PE and QE
            good_hrtf = n_good_pe[n_good_pe['target_HRTF'].isin(n_good_qe['target_HRTF']
                                                                )].target_HRTF.to_numpy()
        # if no HRTFs are classed as good, return a NaN
        else:
            good_hrtf = np.empty(1)
            good_hrtf[:] = np.nan

    # Finding the 'worst' HRTF based on median QE
    if bad_errors_df.empty:
        warnings.warn("No HRTFs classified as bad, selecting the worst one from the good ones.")
        bad_hrtf = good_errors_df.nlargest(1, 'querr_3rd_quartile', keep='all'
                                      ).target_HRTF.to_numpy()
        if len(bad_hrtf) > 1:
            bad_hrtf = good_errors_df[good_errors_df.target_HRTF.isin(bad_hrtf)
            ].nlargest(1, 'rmsP', keep='all').target_HRTF.to_numpy()
    else:
        bad_hrtf = bad_errors_df.nlargest(1, 'querr_median', keep='all'
                                      ).target_HRTF.to_numpy()
        # If multiple HRTFs have the same median QE then check 3rd quartile
        if len(bad_hrtf) > 1:
            bad_hrtf = bad_errors_df[
                bad_errors_df.target_HRTF.isin(bad_hrtf)
            ].nlargest(1, 'querr_3rd_quartile', keep='all').target_HRTF.to_numpy()
            # Finally, check if PEs are different
            if len(bad_hrtf) > 1:
                bad_hrtf = bad_errors_df[
                    bad_errors_df.target_HRTF.isin(bad_hrtf)
                ].nlargest(1, 'rmsP', keep='all').target_HRTF.to_numpy()
    if good_hrtf.size > 1:
        warnings.warn("Multiple HRTFs were classified as best. Returning the first one only")
    if bad_hrtf.size > 1:
        warnings.warn("Multiple HRTFs were classified as worst. Returning the first one only")

    return good_hrtf[0], bad_hrtf[0]


def rms_p(data):
    return np.sqrt(np.sum(np.square(data)) / len(data))





def plot_hrtf_distributions(template: str,
                            target_hrtfs,
                            good_errors_df: pd.DataFrame,
                            best_hrtf: str,
                            bad_errors_df: pd.DataFrame,
                            worst_hrtf: str,
                            dir_errors: pd.DataFrame,
                            save_figure: bool):
    my_pal = {target:
                  my_cmap[2] if target == template else
                  my_cmap[3] if target in
                                good_errors_df.target_HRTF.unique() else
                  my_cmap[1] for target in target_hrtfs}

    pe_order = pd.concat([bad_errors_df.iloc[0:1, 0],  # ind HRTF at the start
                          good_errors_df.sort_values('rmsP').target_HRTF,
                          bad_errors_df.sort_values('rmsP').target_HRTF])

    sns.set_theme(style="whitegrid", font_scale=1.6)
    fg, axes = plt.subplots(2, 2, figsize=[20, 7], layout="constrained",
                            gridspec_kw={'width_ratios': [9, 1], 'wspace': 0.01})

    dir_err_tmp = dir_errors[(dir_errors['template_HRTF'] == template)]

    ax = sns.violinplot(ax=axes[0, 0], data=dir_err_tmp, x="target_HRTF", y="rmsP", scale="width",
                        inner="quartile", split=False, order=pe_order)
    ax.set_ylabel("PE (deg)")
    ax.set_xlabel(None)
    ax.set_xticklabels([i.get_text()[0:5] for i in ax.xaxis.get_ticklabels()])
    ax.set_ylim(bottom=0)

    for ind, violin in enumerate(ax.findobj(PolyCollection)):
        t = pe_order.values[ind]
        violin.set_edgecolor(my_cmap[0] if t == worst_hrtf else
                             my_cmap[4] if t == best_hrtf else
                             'black')
        violin.set_hatch('///' if t == worst_hrtf else
                         '\\\\\\' if t == best_hrtf else
                         None)
        violin.set_facecolor(my_pal[t])

    ax.legend(handles=handles,
              labels=['Limited distributions:', 'Bad', 'Worst', 'Good',
                      'Best', 'Individual'], loc='upper center', bbox_to_anchor=(0.5, 1.15), fontsize=13,
              ncol=7, frameon=False, columnspacing=1.4, handlelength=1.4, handletextpad=0.5)

    qe_order = pd.concat([bad_errors_df.iloc[0:1, 0],  # ind HRTF at the start
                          good_errors_df.sort_values('querr_3rd_quartile').target_HRTF,
                          bad_errors_df.sort_values(['querr_median', 'querr_3rd_quartile']).target_HRTF])
    ax = sns.violinplot(ax=axes[1, 0], data=dir_err_tmp, x="target_HRTF", y="querr", cut=0,
                        scale="width", inner="quartile", order=qe_order, legend=False)
    ax.set_ylabel("QE (%)")
    ax.set_xlabel("Target HRTF")
    ax.set_xticklabels([i.get_text()[0:5] for i in ax.xaxis.get_ticklabels()])
    ax.set_ylim(bottom=0)

    for ind, violin in enumerate(ax.findobj(PolyCollection)):
        t = qe_order.values[ind]
        violin.set_edgecolor(my_cmap[0] if t == worst_hrtf else
                             my_cmap[4] if t == best_hrtf else
                             'black')
        violin.set_hatch('///' if t == worst_hrtf else
                         '\\\\\\' if t == best_hrtf else
                         None)
        violin.set_facecolor(my_pal[t])

    ax.legend([], [], frameon=False)

    sns.set_theme(style="whitegrid", font_scale=0.9)
    dir_err_gb = dir_errors[(dir_errors['template_HRTF'] == template) &
                            ((dir_errors['target_HRTF'] == best_hrtf) |
                             (dir_errors['target_HRTF'] == worst_hrtf))]
    ax = sns.violinplot(ax=axes[0, 1], data=dir_err_gb, x="target_HRTF", y="rmsP", scale="count", inner="quartile",
                        split=False, order=[best_hrtf, worst_hrtf], palette=[my_cmap[-2], my_cmap[1]])
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_xticklabels([i.get_text()[0:5] for i in ax.xaxis.get_ticklabels()])
    ax.tick_params(axis='y', which='major', pad=0)

    for ind, violin in enumerate(ax.findobj(PolyCollection)):
        violin.set_edgecolor(my_cmap[4 - 4 * ind])
        violin.set_hatch('\\\\\\' if ind == 0 else '///')
    ax = sns.violinplot(ax=axes[1, 1], data=dir_err_gb, x="target_HRTF", y="querr", scale="count", inner="quartile",
                        split=False, order=[best_hrtf, worst_hrtf], cut=0, palette=[my_cmap[-2], my_cmap[1]])
    for ind, violin in enumerate(ax.findobj(PolyCollection)):
        violin.set_edgecolor(my_cmap[4 - 4 * ind])
        violin.set_hatch('\\\\\\' if ind == 0 else '///')
    ax.set_ylim(bottom=0)
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    ax.set_xticklabels([i.get_text()[0:5] for i in ax.xaxis.get_ticklabels()])
    ax.tick_params(axis='y', which='major', pad=0)
    if save_figure:
        plt.savefig(Path('selection_figures')/Path(template + '_distribution.pdf'),
                    format='pdf', bbox_inches='tight')
        plt.close(fg)

def selection_matrix(selected_hrtf_list, template_hrtfs, target_hrtfs, good_err_all, bad_err_all):
    zd = np.zeros(shape=(len(target_hrtfs), len(template_hrtfs)))
    sel_m = pd.DataFrame(data=zd, index=target_hrtfs, columns=template_hrtfs)
    for t in template_hrtfs:
        if not selected_hrtf_list[selected_hrtf_list.id == t]['goodHRTF'].isnull().values.any():
            sel_m.loc[good_err_all[good_err_all.template_HRTF == t].target_HRTF, t] = 1
            sel_m.loc[selected_hrtf_list[selected_hrtf_list.id == t]['goodHRTF'], t] = 2

        sel_m.loc[bad_err_all[bad_err_all.template_HRTF == t].target_HRTF, t] = -1
        sel_m.loc[selected_hrtf_list[selected_hrtf_list.id == t]['badHRTF'], t] = -2

    sel_m.columns = [i[0:5] for i in template_hrtfs]
    sel_m.index = [i[0:5] for i in target_hrtfs]
    return sel_m


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'no_plot':
        select_HRTFs(False)
        print('Figures not plotted.')
    else:
        select_HRTFs(True)
        print('Figures plotted.')