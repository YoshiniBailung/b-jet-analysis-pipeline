import os

from src.training import train_model
from src.evaluation import evaluate_validation, evaluate_data, evaluate_performance
from src.plot import plot_score, plot_shap, plot_2d, plot_2d_Db, plot_Db, plot_distributions, plot_correlation
from src.data import load_parquet, cleanup
from src.histograms import save_2d_histograms, save_2d_db_histograms, save_1d_histograms, save_db_histograms
# from src.roofit import load_ali_env, run_roofit_1d, run_roofit_2d


def run_pt_bin(pt_bin, config):

    low, high = pt_bin
    pt_bin_string = f'pt_{low}_{high}'

    train_file = config["paths"]["train_pattern"].format(low=low, high=high)
    data_file = config["paths"]["data_pattern"].format(low=low, high=high)

    os.makedirs(config["paths"]["figure_dir"].format(low=low, high=high), exist_ok = True)
    figure_path = config["paths"]["figure_dir"].format(low=low, high=high)
    
    os.makedirs(config["paths"]["hist_dir"].format(low=low, high=high), exist_ok = True)
    hist_path = config["paths"]["hist_dir"].format(low=low, high=high)


    os.makedirs(config["paths"]["fit_dir"].format(low=low, high=high), exist_ok = True)
    fit_path = config["paths"]["fit_dir"].format(low=low, high=high)
    
    train_df = load_parquet(train_file)
    #train_df = cleanup(train_df)

    plot_distributions(train_df, config, pt_bin, figure_path)
    plot_correlation(train_df, config, pt_bin, figure_path)
    
    model, train_test_data = train_model(train_df, pt_bin_string, config)

    scores = evaluate_validation(model, train_test_data[2], train_test_data[3])

    plot_score(model, train_test_data, pt_bin, figure_path)
    plot_shap(model, train_test_data[2], train_test_data[3], pt_bin, figure_path)

    evaluate_performance(scores, train_test_data[3], pt_bin, figure_path)

    data_df = load_parquet(data_file)
    #data_df = data_df[data_df["nTracks"] >= 3]

    scores_data = evaluate_data(model, data_df, config)

    plot_2d(scores, scores_data, train_test_data[3], pt_bin, figure_path)
    plot_2d_Db(scores, scores_data, train_test_data[3], pt_bin, figure_path)
    plot_Db(scores, scores_data, train_test_data[3], pt_bin, figure_path)

    save_2d_histograms(scores, scores_data, train_test_data[3], hist_path)
    save_2d_db_histograms(scores, scores_data, train_test_data[3], hist_path)
    save_1d_histograms(scores, scores_data, train_test_data[3], hist_path)
    save_db_histograms(scores, scores_data, train_test_data[3], hist_path)
