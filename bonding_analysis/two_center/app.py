import os
import argparse

import dash
from monty.os.path import zpath
from pymatgen.core import Structure
from pymatgen.io.lobster.outputs import (
    Charge,
    Icohplist,
)
import crystal_toolkit.components as ctc
from crystal_toolkit.helpers.layouts import H3, Container
from crystal_toolkit.settings import SETTINGS
from lobsterenv_vis import LobsterEnvComponent


def main(calc_dir: str, read_icobis: bool):
    # Load correct file based on flag
    icohplist_file = zpath(
        os.path.join(
            calc_dir, "ICOBILIST.lobster" if read_icobis else "ICOHPLIST.lobster"
        )
    )

    icohplist_obj = Icohplist(
        filename=icohplist_file, are_cobis=read_icobis, are_coops=False
    )

    charge_obj = Charge(filename=zpath(os.path.join(calc_dir, "CHARGE.lobster")))
    structure_obj = Structure.from_file(zpath(os.path.join(calc_dir, "CONTCAR")))

    lob_env_component = LobsterEnvComponent(
        charge_obj=charge_obj,
        icohplist_obj=icohplist_obj,
        structure_obj=structure_obj,
        disable_callbacks=False,
    )

    type_pop = "ICOHPs" if not icohplist_obj.are_cobis else "ICOBIs"
    app_header = f"{structure_obj.composition.reduced_formula}: {type_pop}"

    layout = Container([H3(app_header), lob_env_component.layout()])
    app = dash.Dash(assets_folder=SETTINGS.ASSETS_PATH, prevent_initial_callbacks=True)
    ctc.register_crystal_toolkit(app, layout=layout)

    app.run(debug=True, port=8051)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize LobsterEnv environments.")
    parser.add_argument("calc_path", type=str, help="Path to LOBSTER calculation directory")
    parser.add_argument(
        "--read-icobis",
        action="store_true",
        help="Read ICOBILIST.lobster instead of default ICOHPLIST.lobster",
    )

    args = parser.parse_args()
    main(args.calc_path, args.read_icobis)
