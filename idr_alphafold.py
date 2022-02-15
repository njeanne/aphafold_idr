#! /usr/bin/env python3

"""
Created on 18 Jan. 2022
"""

import argparse
import csv
import logging
import os
import statistics
import sys

import altair as alt
from altair_saver import save
from Bio.PDB.PDBParser import PDBParser
import pandas as pd

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[args.log_level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def restricted_int(int_to_inspect):
    """Inspect if an int is between 0 and 100

    :param int_to_inspect: the int to inspect
    :type int_to_inspect: int
    :raises ArgumentTypeError: is not between 0 and 100
    :return: the int value if float_to_inspect is between 0 and 100
    :rtype: int
    """
    x = int(int_to_inspect)
    if x < 0 or x > 100:
        raise argparse.ArgumentTypeError("{} not in range [0, 100]".format(x))
    return x


def extract_plddt(input_path):
    """
    Extract the pLDDT (Beta factor column) of the Alphafold PDB file by amino acid.

    :param input_path: the pdb file path
    :type input_path: str
    :return: the residues and pLDDT dataframe.
    :rtype: Pandas.DataFrame
    """
    with open(input_path, "r") as input_pdb:
        plddt_data = {"position": [], "pLDDT": [], "amino-acid": []}
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure("old", input_pdb)
        for model in structure:
            for chain in model:
                for residue in chain.get_residues():
                    plddt_data["position"].append(residue.full_id[3][1])
                    plddt_data["pLDDT"].append(list(residue.get_atoms())[0].get_bfactor())
                    plddt_data["amino-acid"].append(residue.resname)
        return pd.DataFrame(plddt_data)


def get_residue_order_state(plddt_data, threshold, window):
    """
    On window compute the mean of pLDDT values and determine if the region is ordered or disordered from the center of
    the window to each side of the window.

    :param plddt_data: the dataframe of the residues pLDDT values.
    :type plddt_data: Pandas.DataFrame
    :param threshold: the threshold for a disordered region, under or equal.
    :type threshold: float
    :param window: the window size.
    :type window: int
    :return: the updated dataframe of the residues mean pLDDT on the window which center center is the residue position
    and if the residue normalized on the window is ordered or not.
    :rtype: Pandas.DataFrame
    """
    plddt_means = []
    plddt_list = list(plddt["pLDDT"])
    if window % 2 == 0:
        idx_start = int(window / 2) - 1
        left_window = int(window / 2) - 1
        right_window = int(window / 2)
    else:
        idx_start = int(window / 2)
        left_window = int(window / 2)
        right_window = int(window / 2)

    logging.info("window size for disordered regions search: {}".format(window))
    logging.info("threshold for disordered regions <= {}%".format(threshold))
    logging.debug("index start for disordered regions search: {}".format(idx_start))
    logging.debug("window left size for disordered regions search: {}".format(left_window))
    logging.debug("window right size for disordered regions search: {}".format(right_window))

    for idx in range(idx_start):
        mean_plddt = statistics.mean(plddt_list[0:window])
        plddt_means.append(mean_plddt)
        logging.debug(
            "pLDDT mean value computed on the first {} residues (window size), residue {}: {}".format(window, idx + 1,
                                                                                                      mean_plddt))
    for idx in range(idx_start, len(plddt_list) - right_window):
        mean_plddt = statistics.mean(plddt_list[(idx - left_window):(idx + right_window)])
        plddt_means.append(mean_plddt)
        logging.debug("pLDDT mean value computed on the window size ({}), residue {}: {}".format(window, idx + 1,
                                                                                                 mean_plddt))
    for idx in range(len(plddt_list) - right_window, len(plddt_list)):
        mean_plddt = statistics.mean(plddt_list[len(plddt_list)-window:len(plddt_list)])
        plddt_means.append(mean_plddt)
        logging.debug(
            "pLDDT mean value computed on the last {} residues (window size), residue {}: {}".format(window, idx + 1,
                                                                                                     mean_plddt))
    plddt_data["pLDDT window mean"] = plddt_means
    plddt_data["order state"] = ["ordered" if value >= threshold else "disordered" for value in plddt_means]

    return plddt_data


def get_domains(plddt_data, domains):
    """
    Get the domains which each residue belongs to using its position.

    :param plddt_data: the dataframe of the pLDDT.
    :type plddt_data: Pandas.DataFrame
    :param domains: the dataframe of the domains of the protein.
    :type domains: Pandas.DataFrame
    :return: the updated dataframe with the domain which each residue belongs to.
    :rtype: Pandas.DataFrame
    """

    domains_list = []

    for pos in plddt_data["position"]:
        domain = None
        for _, row in domains.iterrows():
            if row["start"] <= pos <= row["end"]:
                domain = row["domain"]
                break
        domains_list.append(domain)

    plddt_data["domains"] = domains_list
    return plddt_data


def get_areas_order_state(plddt_data):
    """
    Get coordinates of the distinct areas depending on their order state.

    :param plddt_data: the dataframe of the pLDDT.
    :type plddt_data: Pandas.DataFrame
    :return: the dataframe of the areas depending on their order state.
    :rtype: Pandas.DataFrame
    """
    order_state_areas = {"state": [], "start": [], "end": [], "color": []}
    state_color = {"ordered": "#1500ff", "disordered": "#ff0000"}
    residue_order_state = None
    for _, row in plddt_data.iterrows():
        if not residue_order_state == row["order state"]:
            if residue_order_state is not None:
                order_state_areas["end"].append(previous_position)
            order_state_areas["state"].append(row["order state"])
            order_state_areas["start"].append(row["position"])
            order_state_areas["color"].append(state_color[row["order state"]])
        residue_order_state = row["order state"]
        previous_position = row["position"]
    # record the last area end position
    order_state_areas["end"].append(row["position"])

    return pd.DataFrame(order_state_areas)


def draw_chart_plddt(plddt_data, threshold, out_dir, prot_id, out_format, domains=None):
    """
    Draw the chart for the pLDDT values.

    :param plddt_data: the dataframe of the pLDDT.
    :type plddt_data: Pandas.DataFrame
    :param threshold: the order / disorder pLDDT threshold.
    :type threshold: int
    :param out_dir: the chart directory output path.
    :type out_dir: str
    :param prot_id: the name of the protein.
    :type prot_id: str
    :param out_format: the format of the output chart file.
    :type out_format: str
    :param domains: the dataframe of the domains of the protein.
    :type domains: Pandas.DataFrame
    :return: the chart directory output path.
    :rtype: str
    """
    out_path = os.path.join(out_dir, "pLDDT_{}.{}".format(prot_id, out_format))

    plddt_chart = alt.Chart(data=plddt_data,
                            title="{}: pLDDT by residue position".format(prot_id)).mark_line(color="#000000",
                                                                                             strokeWidth=0.4).encode(
        x=alt.X("position", title="amino-acids positions"),
        y=alt.Y("pLDDT", title="pLDDT"),
    )

    areas_order_state = get_areas_order_state(plddt_data)
    disorder_chart = alt.Chart(areas_order_state).mark_rect(opacity=0.3).encode(
        x=alt.X("start"),
        x2=alt.X2("end"),
        color=alt.Color("state",
                        scale=alt.Scale(domain=areas_order_state["state"].tolist(),
                                        range=areas_order_state["color"].tolist()))
    )

    # add the order/disorder threshold
    threshold_chart = alt.Chart(pd.DataFrame([{"threshold": threshold}])).mark_rule(color="red").encode(
        y=alt.Y("threshold", title="pLDDT")
    )

    plddt_chart = disorder_chart + threshold_chart + plddt_chart

    if domains is not None:
        domain_chart = alt.Chart(domains).mark_rect().encode(
            x=alt.X("start", title="domains"),
            x2=alt.X2("end"),
            color=alt.Color("domain",
                            scale=alt.Scale(domain=df_domains["domain"].tolist(), range=df_domains["color"].tolist()))
        )
        plddt_chart = alt.vconcat(plddt_chart, domain_chart).resolve_scale(color="independent")

    save(plddt_chart, out_path)
    return out_path


if __name__ == "__main__":
    descr = """
    {} v. {}

    Created by {}.
    Contact: {}
    {}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or
    implied.

    Intrinsically Disordered Regions visualisation on Alphafold outputs, based on the pLDTT (predicted Local Distance 
    Difference Test) values in the .pdb files of the predicted structures.
    See Jumper et al. 2021, Suppl. Methods 1.9.6 and Mariani et al. 2013 Bioinformatics for details.
    
    An optional Comma Separated Values file can be provided with the AA coordinates of the different regions. The file
    must have a formatted header with the following values: \"domain,start,end,color\".
    - domain:   the domain name
    - start:    the domain AA start (1-index)
    - end:      the domain AA end (1-index)
    - color:    the color used to display the domain (hexadecimal format)
    """.format(os.path.basename(__file__), __version__, __author__, __email__, __copyright__)
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="path to the output directory.")
    parser.add_argument("-f", "--format", required=False, type=str, default="svg",
                        choices=["svg", "html"], help="the output format of the plot, default is SVG.")
    parser.add_argument("-d", "--domains", required=False, type=str,
                        help="CSV (comma separated) file defining the domains name, the start position, the stop "
                             "position and the color representation in hexadecimal format.")
    parser.add_argument("-w", "--window-size", required=False, type=int, default=11,
                        help="the window size to determine if the mean of the pLDDT values in the window are ordered "
                             "or disordered, default is 11.")
    parser.add_argument("-t", "--threshold", required=False, type=restricted_int, default=50,
                        help="the threshold percentage of disorder for plDDT. If the pLDDT is under or equal to this "
                             "threshold the pLDDT value is set as disordered. Default value is 50%%.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help=("the path for the log file. If this option is  skipped, the log file is created in the "
                              "output directory."))
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If this option is skipped, the log level is INFO.")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument("input", type=str, help="path to the Alphafold prediction .pdb file.")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, "{}.log".format(os.path.splitext(os.path.basename(__file__))[0]))
    create_log(log_path, args.log_level)

    logging.info("version: {}".format(__version__))
    logging.info("CMD: {}".format(" ".join(sys.argv)))

    alphafold_prediction_id = os.path.splitext(os.path.basename(args.input))[0]
    plddt = extract_plddt(args.input)
    plddt = get_residue_order_state(plddt, float(args.threshold), args.window_size)

    if args.domains:
        df_domains = pd.read_csv(args.domains, sep=",", header=0)
        plddt = get_domains(plddt, df_domains)
        path_chart = draw_chart_plddt(plddt, args.threshold, args.out, alphafold_prediction_id, args.format, df_domains)
    else:
        path_chart = draw_chart_plddt(plddt, args.threshold, args.out, alphafold_prediction_id, args.format)

    logging.info("pLDDT chart for {} created: {}".format(alphafold_prediction_id, os.path.abspath(path_chart)))
    path_data = os.path.join(args.out, "pLDDT_{}.csv".format(alphafold_prediction_id))
    plddt.to_csv(path_data, index=False)
    logging.info("pLDDT data for {} created: {}".format(alphafold_prediction_id, os.path.abspath(path_data)))
