import subprocess

def load_ali_env():

    cmd = 'alienv enter O2Physics/latest-master-o2'

    subprocess.run(cmd, shell = True, check = True)

def run_roofit_1d(fit_path):

    macro = config["roofit"]["macro1"]
    cmd = f'cp {macro} fit_path && cd fit_path && root -l -b -q {macro}("data_1d.root", "hist", "template_1d_b.root", "hist_b", "template_1d_c.root", "hist_c", "template_1d_lf.root", "hist_lf")'

    subprocess.run(cmd, shell=True, check=True)

def run_roofit_2d(fit_path):

    macro = config["roofit"]["macro2"]

    cmd = f'cp {macro} fit_path && cd fit_path && root -l -b -q {macro}("data_2d.root", "hist", "template_2d_b.root", "hist_b", "template_2d_c.root", "hist_c", "template_2d_lf.root", "hist_lf")'

    subprocess.run(cmd, shell=True, check=True)