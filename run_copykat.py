# running copycat
# some info --> copykat only uses hg38 (called hg20 here)

# imports
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
import rpy2.robjects as robjects
from pathlib import Path
import shutil
from pyomics.utils import benchmark_method
import itertools
# ----------------------------------------------------------------------------------------------------------------------

def grid_by_dict(pars_dict: dict) -> list:
    keys=pars_dict.keys()
    combinations = itertools.product(*pars_dict.values())
    list_of_kwargs = [dict(zip(keys, cc)) for cc in combinations]
    return list_of_kwargs

def val_build_project() -> (Path, Path):
    cwd_path = Path.cwd()
    print(f"Current working directory of running script {Path(__file__).name}: {cwd_path}")
    path_out = cwd_path / "app" / "out"
    path_in = cwd_path / "data_input"

    if not path_in.exists():
        raise ValueError(f"Data dir '{str(path_in)}' does not exist!")

    if not path_out.exists():
        path_out.mkdir(parents=True, exist_ok=True)
        print(f"Data out-dir '{str(path_out)}' has been created...")

    return path_in, path_out


def get_hg_38_file_paths(target_path: Path) -> list:
    return [p for p in target_path.rglob("*__hg_38__RCM.csv") if p.is_file()]


def get_hg_38_desc_paths(target_path: Path) -> dict:
    return {p.stem: p for p in target_path.rglob("*__hg_38.txt")}


def run_copykat(path_target: Path,
                path_out_data: Path,
                n_cores: int = 10,
                n_genes_chr: int = 5,
                window_size: int = 25,
                low_dr: float = 0.05,
                up_dr: float = 0.1,
                ks_cut: float = 0.1,
                cell_pre_label: bool = False):
    list_paths_target_csvs = get_hg_38_file_paths(path_target)
    if cell_pre_label:
        dict_paths_target_txts = get_hg_38_desc_paths(path_target)
        list_paths_target_csvs = [p for p in list_paths_target_csvs if p.stem.split("__RCM")[0] in dict_paths_target_txts.keys()]
    for p in list_paths_target_csvs:
        name_tag = f"{p.stem}__n,{n_cores};n,{n_genes_chr};w,{window_size};l,{low_dr};u,{up_dr};k,{ks_cut};c,{cell_pre_label}__copykat_"
        path_out_target = path_out_data / f"out__{name_tag}"
        path_out_target.mkdir(parents=True, exist_ok=True)

        norm_cell_vector = ""
        if cell_pre_label:
            path_txt = dict_paths_target_txts[p.stem.split("__RCM")[0]]   # normal cells split by \n
            with open(path_txt, "r") as f:
                list_norm_cells = list(map(lambda x: x.replace("\n", ""), f.readlines()))
                norm_cell_vector = robjects.vectors.StrVector(list_norm_cells)

        ################################################################################################################
        # run copykat R-script
        @benchmark_method(path_out_target)
        def run_rscript(p,
                        name_tag,
                        n_cores,
                        n_genes_chr,
                        window_size,
                        ks_cut,
                        low_dr,
                        up_dr,
                        norm_cell_vector):
            r = robjects.r
            r.source("c_copykatR.R")
            r.r_run_copykat(str(p), name_tag, n_cores, n_genes_chr, window_size, ks_cut, low_dr, up_dr, norm_cell_vector)

        run_rscript(p,
                    name_tag,
                    n_cores,
                    n_genes_chr,
                    window_size,
                    ks_cut,
                    low_dr,
                    up_dr,
                    norm_cell_vector)
        ################################################################################################################

        # reformat *__copykat_CNA_raw_results_gene_by_cell.txt to final GBC-format
        path_pre_gbc = [p for p in Path.cwd().glob("*__copykat_CNA_raw_results_gene_by_cell.txt")][0]
        df_gbc_pre = pd.read_csv(path_pre_gbc, sep="\t")
        df_gbc_export = df_gbc_pre.drop(["abspos", "band", "ensembl_gene_id", "hgnc_symbol"], axis=1).rename(
                 {"chromosome_name": "CHR", "start_position": "START", "end_position": "END"}, axis=1).set_index("CHR")
        df_gbc_export.to_csv(Path.cwd() / f"{name_tag}_copykat__GBC.csv")

        # export everything to /app/out
        list_all_nametag_items = [p for p in Path.cwd().glob(f"*{name_tag}*")]
        for items in list_all_nametag_items:
            shutil.move(items, path_out_target / items.name)


if __name__ == "__main__":

    # matrix of possible copykat hyperparameter kwargs
    kwargs_gridsearch = {
        "n_cores": [5, 10, 20],
        "n_genes_chr": [1, 5, 10, 100],
        "window_size": [10, 25, 50, 100, 200, 500],
        "low_dr": [0, 0.05, 0.1, 0.2],
        "up_dr": [0, 0.1, 0.2],
        "ks_cut": [0, 0.1, 0.2],
        "cell_pre_label": [True, False]
    }

    path_in, path_out = val_build_project()
    run_copykat(path_in, path_out, n_cores=1, cell_pre_label=False)  # standard parameters; one core
    run_copykat(path_in, path_out, n_cores=1, cell_pre_label=True)
    list_kwargs = grid_by_dict(kwargs_gridsearch)
    for kwarg_opt in list_kwargs:
        print(f"CopyKat running with hyperparameters: {kwarg_opt}")
        run_copykat(path_in, path_out, **kwarg_opt)
