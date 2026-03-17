import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

def evaluate_validation(model, X_test, y_test):

    scores = model.predict(X_test, output_margin = False)

    roc_auc = roc_auc_score(y_test, scores, multi_class='ovo')

    print(f"ROC AUC Score: {roc_auc}")

    return scores

def evaluate_data(model, df, config):

    features = config["features"]

    X = df[features]

    scores = model.predict(X, output_margin=False)

    return scores

def evaluate_performance(scores, y_test, pt_bin, path):

    low, high = pt_bin

    abundances = [1,2,20]
    
    if (low == 5):
        abundances = [1,3,43]
    elif (low >= 10 & high <= 20):
        abundances = [1,2,25]
    else:
        abundances = [1,2,20]

    test_scores_lf = scores[y_test == 0]
    test_scores_c = scores[y_test == 1]
    test_scores_b = scores[y_test == 2]

    #true b-jet discriminator
    b_score_b_c = test_scores_b[:,2] / (test_scores_b[:, 1] + test_scores_b[:, 2]) 
    b_score_bc_usdg = (test_scores_b[:,2] + test_scores_b[:,1])/ (test_scores_b[:, 0] + test_scores_b[:,1] + test_scores_b[:,2]) 

    #true-c-jet discriminator
    c_score_b_c = test_scores_c[:,2] / (test_scores_c[:, 1] + test_scores_c[:, 2])
    c_score_bc_usdg = (test_scores_c[:,2] + test_scores_c[:,1])/ (test_scores_c[:, 0] + test_scores_c[:,1] + test_scores_c[:,2])

    #true-lf-jet discriminator
    lf_score_b_c = test_scores_lf[:,2] / (test_scores_lf[:, 1] + test_scores_lf[:, 2])
    lf_score_bc_usdg = (test_scores_lf[:,2] + test_scores_lf[:,1])/ (test_scores_lf[:, 0] + test_scores_lf[:,1] + test_scores_lf[:,2]) 

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

    w_b = len(test_scores_b)*abundances[0]
    w_c = len(test_scores_c)*abundances[1]
    w_lf = len(test_scores_lf)*abundances[2]

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

    X, Y = np.meshgrid(xedges, yedges)
    Z = significance

    fig, ax = plt.subplots()

    cont = ax.contourf(X, Y, Z.T)

    ax.set_xlabel('BDT (bc|usdg)', fontsize=20)
    ax.set_ylabel('BDT (b|c)', fontsize=18)

    ax.tick_params(axis='both', which='major', labelsize=18)

    fig.colorbar(cont, ax=ax)

    # legend trick for efficiency/purity
    ax.plot([], [], label=f'Efficiency = {eff_at_max_sig:.4f} | Purity = {purity_at_max_sig:.4f}', color='white')

    ax.axhline(y=y_cut,color='red',label=f'b vs c cut = {y_cut:.4f}',linewidth=3)

    ax.axvline(x=x_cut,color='blue',label=f'bc vs usdg cut = {x_cut:.4f}',linewidth=3)

    ax.legend(fontsize=12,title=f'Significance distribution\n[{low}-{high}] GeV/c | MC LHC25a2b',framealpha=1,title_fontsize=12)
    fig.savefig(f"{path}/significance_2D_distrib_{low}_{high}.pdf",bbox_inches='tight')
    plt.close(fig)
