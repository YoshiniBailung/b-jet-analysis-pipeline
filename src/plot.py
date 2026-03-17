
import numpy as np
from hipe4ml import plot_utils
import matplotlib.pyplot as plt

def plot_correlation(df, config, pt_bin, path):
    
    low, high = pt_bin

    features = config["features"]
    leg_labels = ["lf", "c", "b"]
    label = config["label"]
    df_lf = df[df[label] == 0]
    df_c = df[df[label] == 1]
    df_b = df[df[label] == 2]

    df_all = [df_lf, df_c, df_b]

    corr_fig = plot_utils.plot_corr(df_all, features, leg_labels)

    for fig, label in zip(corr_fig, leg_labels):

        ax = fig.axes[0]

        ax.set_title(f"{label} output | [{low}-{high}] GeV/c", fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        fig.savefig(f"{path}/correlation_{label}_{low}-{high}.pdf", bbox_inches="tight")

        plt.close(fig)

def plot_distributions(df, config, pt_bin, path):
    low, high = pt_bin

    features = config["features"]
    leg_labels = ["lf", "c", "b"]
    label = config["label"]
    df_lf = df[df[label] == 0]
    df_c = df[df[label] == 1]
    df_b = df[df[label] == 2]

    df_all = [df_lf, df_c, df_b]

    dist_plot = plot_utils.plot_distr(df_all, features, bins=100, labels=leg_labels, log=True, density=True, figsize=(18,12), alpha=0.3, grid=False)
    
    fig = plt.gcf()

    fig.subplots_adjust(left=0.06,bottom=0.06,right=0.99,top=0.96,hspace=0.55,wspace=0.55)

    fig.savefig(f"{path}/Feature_distributions.pdf",bbox_inches="tight")

    plt.close(fig)

def plot_score(model, train_test_data, pt_bin, path):

    low, high = pt_bin
    leg_labels = ["lf", "c", "b"]

    ml_out_fig = plot_utils.plot_output_train_test(model, train_test_data, bins=100, 
                                               output_margin=False, labels=leg_labels,
                                               logscale=True, density=True)
    labels = ["lf", "c", "b"]

    for fig, label in zip(ml_out_fig, labels):

        ax = fig.axes[0]

        ax.set_title(f"{label} output | [{low}-{high}] GeV/c", fontsize=30)
        ax.tick_params(axis='both', which='major', labelsize=25)
        ax.set_xlim(-0.05,1.05)

        ax.set_xlabel(f"BDT output for {label}", fontsize=30)
        ax.set_ylabel("Counts (arb. units)", fontsize=30)

        fig.savefig(f"{path}/bdt_output_{label}_{low}-{high}.pdf", bbox_inches="tight")

        plt.close(fig)

def plot_shap(model, X_test, y_test, pt_bin, path):

    low, high = pt_bin
    leg_labels = ["lf", "c", "b"]

    shap_figs = plot_utils.plot_feature_imp(X_test, y_test, model, labels=leg_labels)
    fig = shap_figs[0]
    ax = fig.axes[0]

    ax.grid()
    ax.set_title(f"Feature importance | [{low}-{high}] GeV/c", fontsize=45)
    ax.legend(fontsize=40, loc='lower right')
    ax.set_xlim(0,1.4)
    ax.tick_params(axis='both', which='major', labelsize=25)

    ax.set_ylabel('', fontsize=200)
    ax.set_xlabel('mean (|SHAP value|)\n(average impact on model output magnitude)',fontsize=35)

    fig.savefig(f"{path}/Feature_importance_{low}-{high}.pdf",bbox_inches='tight')

    for fig in shap_figs:
        plt.close(fig)

def plot_2d(scores, scores_data, y_test, pt_bin, path):

    low, high = pt_bin

    test_scores_lf = scores[y_test == 0]
    test_scores_c = scores[y_test == 1]
    test_scores_b = scores[y_test == 2]

    #true b-jet discriminator
    b_score_b_c = test_scores_b[:,2] / (test_scores_b[:, 1] + test_scores_b[:, 2]) 
    b_score_bc_usdg = (test_scores_b[:,2] + test_scores_b[:, 1])/ (test_scores_b[:, 0] + test_scores_b[:, 1] + test_scores_b[:,2]) 

    #true-c-jet discriminator
    c_score_b_c = test_scores_c[:,2] / (test_scores_c[:, 1] + test_scores_c[:, 2])
    c_score_bc_usdg = (test_scores_c[:,2] + test_scores_c[:, 1])/ (test_scores_c[:, 0] + test_scores_c[:, 1] + test_scores_c[:,2])

    #true-lf-jet discriminator
    lf_score_b_c = test_scores_lf[:,2] / (test_scores_lf[:, 1] + test_scores_lf[:, 2])
    lf_score_bc_usdg = (test_scores_lf[:,2] + test_scores_lf[:, 1])/ (test_scores_lf[:, 0] + test_scores_lf[:, 1] + test_scores_lf[:,2]) 

    #data-jet discriminator
    data_score_b_c = scores_data[:, 2] / (scores_data[:, 1] + scores_data[:, 2]) 
    data_score_bc_usdg = (scores_data[:, 2] + scores_data[:, 1])/ (scores_data[:, 0] + scores_data[:, 1] + scores_data[:, 2])
    
    fig, ax = plt.subplots()

    ax.hist2d(b_score_bc_usdg, b_score_b_c, bins = 25, cmap = 'Reds')

    ax.set_xlabel('BDT(bc|usdg)', fontsize=20)
    ax.set_ylabel('BDT(b|c)', fontsize=18)
    ax.set_title('b-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(f"{path}/2D-discriminant-bjets-{low}-{high}.pdf",bbox_inches = 'tight')

    fig, ax = plt.subplots()

    ax.hist2d(c_score_bc_usdg, c_score_b_c, bins = 25, cmap = 'Reds')

    ax.set_xlabel('BDT(bc|usdg)', fontsize=20)
    ax.set_ylabel('BDT(b|c)', fontsize=18)
    ax.set_title('c-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(f"{path}/2D-discriminant-cjets-{low}-{high}.pdf",bbox_inches = 'tight')

    fig, ax = plt.subplots()

    ax.hist2d(lf_score_bc_usdg, lf_score_b_c, bins = 25, cmap = 'Reds')

    ax.set_xlabel('BDT(bc|usdg)', fontsize=20)
    ax.set_ylabel('BDT(b|c)', fontsize=18)
    ax.set_title('lf-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(f"{path}/2D-discriminant-lfjets-{low}-{high}.pdf",bbox_inches = 'tight')

    fig, ax = plt.subplots()

    ax.hist2d(data_score_bc_usdg, data_score_b_c, bins = 25, cmap = 'Reds')

    ax.set_xlabel('BDT(bc|usdg)', fontsize=20)
    ax.set_ylabel('BDT(b|c)', fontsize=18)
    ax.set_title('Data-jets 22o', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(f"{path}/2D-discriminant-datajets-{low}-{high}.pdf",bbox_inches = 'tight')

    fig, ax = plt.subplots()

    ax.hist(b_score_bc_usdg, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(bc|usdg)')
    ax.hist(b_score_b_c, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(b|c)')


    ax.set_xlabel('BDT(X|Y)', fontsize=20)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_title('b-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(prop={"size": 15}, loc="upper right")
    plt.savefig(f"{path}/1D-discriminant-bjets-{low}-{high}.pdf",bbox_inches = 'tight')


    fig, ax = plt.subplots()

    ax.hist(c_score_bc_usdg, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(bc|usdg)')
    ax.hist(c_score_b_c, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(b|c)')


    ax.set_xlabel('BDT(X|Y)', fontsize=20)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_title('c-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(prop={"size": 15}, loc="upper right")
    plt.savefig(f"{path}/1D-discriminant-cjets-{low}-{high}.pdf",bbox_inches = 'tight')



    fig, ax = plt.subplots()

    ax.hist(lf_score_bc_usdg, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(bc|usdg)')
    ax.hist(lf_score_b_c, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(b|c)')


    ax.set_xlabel('BDT(X|Y)', fontsize=20)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_title('lf-jets MC 25a2b', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(prop={"size": 15}, loc="upper right")
    plt.savefig(f"{path}/1D-discriminant-lfjets-{low}-{high}.pdf",bbox_inches = 'tight')


    fig, ax = plt.subplots()

    ax.hist(data_score_bc_usdg, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(bc|usdg)')
    ax.hist(data_score_b_c, bins = 100, log = True, density = True, alpha = 0.4, label = 'BDT(b|c)')


    ax.set_xlabel('BDT(X|Y)', fontsize=20)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_title('data-jets 22o', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(prop={"size": 15}, loc="upper right")
    plt.savefig(f"{path}/1D-discriminant-datajets-{low}-{high}.pdf",bbox_inches = 'tight')


def plot_Db(scores, scores_data, y_test, pt_bin, path):

    low, high = pt_bin

    fc = 0.018
    test_scores_lf = scores[y_test == 0]
    test_scores_c = scores[y_test == 1]
    test_scores_b = scores[y_test == 2]

    main_term_b = test_scores_b[:,2]/((1-fc)*test_scores_b[:, 0] + fc*test_scores_b[:, 1])
    Db_b = np.log(main_term_b)

    main_term_c = test_scores_c[:,2]/((1-fc)*test_scores_c[:, 0] + fc*test_scores_c[:, 1])
    Db_c = np.log(main_term_c)

    main_term_lf = test_scores_lf[:,2]/((1-fc)*test_scores_lf[:, 0] + fc*test_scores_lf[:, 1])
    Db_lf = np.log(main_term_lf)

    main_term_data = scores_data[:,2]/((1-fc)*scores_data[:, 0] + fc*scores_data[:, 1])
    Db_data = np.log(main_term_data)

    fig, ax = plt.subplots()

    ax.hist(Db_lf, bins = 100, log = True, density = True, alpha = 0.4, label = 'lf-jets')
    ax.hist(Db_c, bins = 100, log = True, density = True, alpha = 0.4, label = 'c-jets')
    ax.hist(Db_b, bins = 100, log = True, density = True, alpha = 0.4, label = 'b-jets')
    ax.hist(Db_data, bins = 100, log = True, density = True, alpha = 0.4, label = 'data-jets')
    ax.set_xlabel('Db', fontsize=20)
    ax.set_ylabel('Counts', fontsize=18)
    ax.set_title('Db distribution', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.legend(prop={"size": 15}, loc="upper right")
    plt.savefig(f"{path}/Db-discriminant-{low}-{high}.pdf",bbox_inches = 'tight')

