import numpy as numpy
import uproot
import boost_histogram as bh

bins = 100
xmin = 0.
xmax = 1.
ymin = 0.
ymax = 1.
dbbins = 200
dbmin = -10.
dbmax = 10.

def save_2d_histograms(scores, data_scores, y_test, path):
    
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

    #data discriminators
    data_score_b_c = data_scores[:, 2] / (data_scores[:, 1] + data_scores[:, 2]) 
    data_score_bc_usdg = (data_scores[:, 2] + data_scores[:, 1])/ (data_scores[:, 0] + data_scores[:, 1] + data_scores[:, 2]) 
    
    # data_score_b_c += 0.0075
    # data_score_bc_usdg += 0.015

    print(len(b_score_b_c), len(b_score_b_c) + len(c_score_b_c) + len(c_score_b_c))

    hdata = bh.Histogram(bh.axis.Regular(bins, xmin, xmax), bh.axis.Regular(bins, ymin, ymax))
    hdata.fill(data_score_bc_usdg, data_score_b_c)

    root_file = uproot.recreate(f"{path}/data_2d.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist"] = hdata.to_numpy()

    htempb = bh.Histogram(bh.axis.Regular(bins, xmin, xmax), bh.axis.Regular(bins, ymin, ymax))
    htempb.fill(b_score_bc_usdg, b_score_b_c)

    root_file = uproot.recreate(f"{path}/template_2d_b.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_b"] = htempb.to_numpy()

    htempc = bh.Histogram(bh.axis.Regular(bins, xmin, xmax), bh.axis.Regular(bins, ymin, ymax))
    htempc.fill(c_score_bc_usdg, c_score_b_c)

    root_file = uproot.recreate(f"{path}/template_2d_c.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_c"] = htempc.to_numpy()

    htemplf = bh.Histogram(bh.axis.Regular(bins, xmin, xmax), bh.axis.Regular(bins, ymin, ymax))
    htemplf.fill(lf_score_bc_usdg, lf_score_b_c)

    root_file = uproot.recreate(f"{path}/template_2d_lf.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_lf"] = htemplf.to_numpy()

def save_1d_histograms(scores, data_scores, y_test, path):

    data_score_b = data_scores[:, 2]

    scores_b = scores[:, 2]

    scores_b_lf = scores_b[y_test == 0]
    scores_b_c = scores_b[y_test == 1]
    scores_b_b = scores_b[y_test == 2]

    print(len(scores_b_b), len(scores_b_b) + len(scores_b_c) + len(scores_b_lf))

    hdata = bh.Histogram(bh.axis.Regular(bins, xmin, xmax))
    hdata.fill(data_score_b)

    root_file = uproot.recreate(f"{path}/data_1d.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist"] = hdata.to_numpy()

    htempb = bh.Histogram(bh.axis.Regular(bins, xmin, xmax))
    htempb.fill(scores_b_b)

    root_file = uproot.recreate(f"{path}/template_1d_b.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_b"] = htempb.to_numpy()


    htempc = bh.Histogram(bh.axis.Regular(bins, xmin, xmax))
    htempc.fill(scores_b_c)

    root_file = uproot.recreate(f"{path}/template_1d_c.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_c"] = htempc.to_numpy()

    htemplf = bh.Histogram(bh.axis.Regular(bins, xmin, xmax))
    htemplf.fill(scores_b_lf)

    root_file = uproot.recreate(f"{path}/template_1d_lf.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_lf"] = htemplf.to_numpy()

def save_db_histograms(scores, data_scores, y_test, path):

    fc = 0.018
    test_scores_lf = scores[y_test == 0]
    test_scores_c = scores[y_test == 1]
    test_scores_b = scores[y_test == 2]

    main_term_b = test_scores_b[:,2]/((1-fc)*test_scores_b[:, 0] + fc*test_scores_b[:, 1])
    Db_b = numpy.log(main_term_b)

    main_term_c = test_scores_c[:,2]/((1-fc)*test_scores_c[:, 0] + fc*test_scores_c[:, 1])
    Db_c = numpy.log(main_term_c)

    main_term_lf = test_scores_lf[:,2]/((1-fc)*test_scores_lf[:, 0] + fc*test_scores_lf[:, 1])
    Db_lf = numpy.log(main_term_lf)

    main_term_data = data_scores[:,2]/((1-fc)*data_scores[:, 0] + fc*data_scores[:, 1])
    Db_data = numpy.log(main_term_data)

    hdata = bh.Histogram(bh.axis.Regular(dbbins, dbmin, dbmax))
    hdata.fill(Db_data)

    root_file = uproot.recreate(f"{path}/data_Db.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist"] = hdata.to_numpy()

    htempb = bh.Histogram(bh.axis.Regular(dbbins, dbmin, dbmax))
    htempb.fill(Db_b)

    root_file = uproot.recreate(f"{path}/template_Db_b.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_b"] = htempb.to_numpy()


    htempc = bh.Histogram(bh.axis.Regular(dbbins, dbmin, dbmax))
    htempc.fill(Db_c)

    root_file = uproot.recreate(f"{path}/template_Db_c.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_c"] = htempc.to_numpy()

    htemplf = bh.Histogram(bh.axis.Regular(dbbins, dbmin, dbmax))
    htemplf.fill(Db_lf)

    root_file = uproot.recreate(f"{path}/template_Db_lf.root")
    root_file["words"] = "See what is in the ROOT File!"
    root_file["hist_lf"] = htemplf.to_numpy()
