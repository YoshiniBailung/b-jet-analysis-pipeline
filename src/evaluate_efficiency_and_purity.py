import numpy as np
import matplotlib.pyplot as plt

def evaluate_efficiency_and_purity(scores, y_test, config):

    test_scores_lf = scores[y_test == 0]
    test_scores_c = scores[y_test == 1]
    test_scores_b = scores[y_test == 2]

    pt_bin = config["pt_bin"]

    low, high = pt_bin

    if (low == 5):
        abundances = [1,3,43]
    elif (low >= 10 && high <= 20):
        abundances = [1,2,25]
    else:
        abundances = [1,2,20]

    #true b-jet discriminator
    b_score_b_c = test_scores_b[:,2] / (test_scores_b[:, 1] + test_scores_b[:, 2]) 
    b_score_bc_usdg = (test_scores_b[:,2])/ (test_scores_b[:, 0] + test_scores_b[:,2]) 

    #true-c-jet discriminator
    c_score_b_c = test_scores_c[:,2] / (test_scores_c[:, 1] + test_scores_c[:, 2])
    c_score_bc_usdg = (test_scores_c[:,2])/ (test_scores_c[:, 0] + test_scores_c[:,2])

    #true-lf-jet discriminator
    lf_score_b_c = test_scores_lf[:,2] / (test_scores_lf[:, 1] + test_scores_lf[:, 2])
    lf_score_bc_usdg = (test_scores_lf[:,2])/ (test_scores_lf[:, 0] + test_scores_lf[:,2]) 

    #2D Efficiency calculation

    bjet_2dhist, xedges, yedges = np.histogram2d(b_score_bc_usdg, b_score_b_c, bins = 200, range = [[0, 1], [0, 1]])
    cjet_2dhist, _ , _ = np.histogram2d(c_score_bc_usdg, c_score_b_c, bins = 200, range = [[0, 1], [0, 1]])
    lfjet_2dhist, _, _ = np.histogram2d(lf_score_bc_usdg, lf_score_b_c, bins = 200, range = [[0, 1], [0, 1]])

    n_pts = 201
    x_threshold = np.linspace(0,1,n_pts)
    y_threshold = np.linspace(0,1,n_pts)

    delta = 1/200.
    matched_efficiency = 0.44

    significance = np.zeros((n_pts, n_pts))
    efficiency = np.zeros((n_pts, n_pts))
    efferror = np.zeros((n_pts, n_pts))
    purity = np.zeros((n_pts, n_pts))
    purerror = np.zeros((n_pts, n_pts))

    weight_lf = 1./0.15
    r_lf = 0.15

    w_b = len*(test_scores_b)*abundances[2]
    w_c = len*(test_scores_c)*abundances[1]
    w_lf = len*(test_scores_lf)*abundances[0]

    for i in range(significance.shape[0]):

        x_bin = np.searchsorted(xedges, x_threshold[i])

        for j in range(significance.shape[1]):

        
            y_bin = np.searchsorted(yedges, y_threshold[j])

            b_sel = bjet_2dhist[x_bin:, y_bin:]

            eff_numerator = np.sum(b_sel)
            eff_denominator = np.sum(bjet_2dhist)#len(test_scores_b[:,2])
            efficiency[i,j] = eff_numerator / eff_denominator
            #efferror[i, j] = np.sqrt((efficiency[i,j]*(1 - efficiency[i,j]))/eff_denominator)

            c_sel = cjet_2dhist[x_bin:, y_bin:]
            lf_sel = lfjet_2dhist[x_bin:, y_bin:]

            purity_denominator = w_b*np.sum(b_sel) + w_c*np.sum(c_sel) + w_lf*np.sum(lf_sel)
            purity_numerator = w_b*np.sum(b_sel)

            if (np.isnan((purity_numerator/purity_denominator))):
                purity[i,j] = 0
                purerror[i,j] = 0
            else:
                purity[i,j] = purity_numerator/purity_denominator
                purerror[i,j] = np.sqrt(np.power((purity_denominator - purity_numerator)/(purity_denominator *purity_denominator),2)*purity_numerator + np.power(purity_numerator/(purity_denominator * purity_denominator),2)*(np.sum(c_sel) + np.sum(lf_sel)*weight_lf*weight_lf) )

            if (efficiency[i,j] < matched_efficiency + delta and efficiency[i,j] > matched_efficiency - delta):
                print(efficiency[i,j], purity[i,j])

            sig_denominator = np.sqrt(w_b*np.sum(b_sel) + w_c*np.sum(c_sel) + w_lf*np.sum(lf_sel))
            sig_numerator = w_b*np.sum(b_sel)

            if (np.isnan(sig_numerator/sig_denominator)):
                significance[i,j] = 0
            else:
                significance[i,j] = sig_numerator/sig_denominator


    max_significance = significance.max()
    index_of_max = np.where(significance == max_significance)
    print(max_significance, index_of_max)
    x_cut = xedges[index_of_max[0][0]]
    y_cut = yedges[index_of_max[1][0]]

    x_cut_bin = np.searchsorted(xedges, x_cut)
    y_cut_bin = np.searchsorted(yedges, y_cut)

    b_sel_at_max_sig = bjet_2dhist[x_cut_bin:, y_cut_bin:]

    eff_at_max_sig = np.sum(b_sel_at_max_sig)/ np.sum(bjet_2dhist)

    c_sel_at_max_sig = cjet_2dhist[x_cut_bin:, y_cut_bin:]
    lf_sel_at_max_sig = lfjet_2dhist[x_cut_bin:, y_cut_bin:]

    purity_at_max_sig = w_b* np.sum(b_sel_at_max_sig)/ (w_b*np.sum(b_sel_at_max_sig) + w_c*np.sum(c_sel_at_max_sig) + w_lf*np.sum(lf_sel_at_max_sig))

    print(x_cut, y_cut, max_significance, eff_at_max_sig, purity_at_max_sig)
    print(efferror[x_cut_bin, y_cut_bin], purerror[x_cut_bin, y_cut_bin])

    plt.figure(1)

    X, Y = np.meshgrid(xedges, yedges)
    Z = significance
    plt.contourf(X, Y, Z.T)
    plt.xlabel('BDT (bc|usdg)', fontsize=20)
    plt.ylabel('BDT (b|c)', fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.colorbar()

    fake_input = np.array([])
    xlabels_input = []
    plt.plot(xlabels_input,fake_input, label='Efficiency = %.4f | Purity = %.4f' % (eff_at_max_sig,purity_at_max_sig), color='white')

    plt.axhline(y = y_cut, color = 'red', label = 'b vs c cut = %.4f' % y_cut,linewidth=3)
    plt.axvline(x = x_cut, color = 'blue', label = 'bc vs usdg cut = %.4f' % x_cut, linewidth=3)

    plt.legend(fontsize=12, title='Significance distribution \n[20-40] GeV/c | MC LHC25a2b',framealpha=1,title_fontsize=12)
    plt.savefig(f{"significance_2D_distrib_{pt_bin}.png"},bbox_inches = 'tight')
