from __future__ import annotations

import os
import pandas as pd
from dash.dependencies import Component, Input, Output
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.core import Molecule, Structure
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.io.lobster.outputs import (
    Charge,
    Icohplist,
)
from pymatgen.util.string import unicodeify_species

from lobsterpy.cohp.describe import Description

from crystal_toolkit.components.structure import StructureMoleculeComponent
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import (
    H4,
    Column,
    Columns,
    get_table,
    html,
)


class LobsterEnvComponent(MPComponent):
    def __init__(
        self,
        charge_obj: Charge | None = None,
        icohplist_obj: Icohplist | None = None,
        structure_obj: Structure | None = None,
        id: str | None = None,
        **kwargs,
    ) -> None:
        super().__init__(
            id=id,
            default_data={
                "charge_obj": charge_obj,
                "icohplist_obj": icohplist_obj,
                "structure_obj": structure_obj,
            },
            **kwargs,
        )

    @property
    def _sub_layouts(self) -> dict[str, Component]:

        data = LobsterEnvComponent._get_all_inputs(self.initial_data["default"])

        analysis_options = [
            {"label": "all", "value": "all"},
            {"label": "cation-anion", "value": "cation-anion"},
        ]

        state = {"analysis-mode": "all"}

        analysis_mode = html.Div(
            [
                self.get_choice_input(
                    kwarg_label="analysis-mode",
                    state=state,
                    label="LobsterEnv analysis mode",
                    help_str="Analysis mode to choose from",
                    options=analysis_options,
                )
            ],
            style={"width": "200px"},
            id=self.id("options-container"),
        )

        # LobsterEnv local environments
        local_envs = html.Div(
            children=[
                LobsterEnvComponent.get_lobster_local_envs(
                    charge_obj=data.get("charge_obj"),
                    icohplist_obj=data.get("icohplist_obj"),
                    structure_obj=data.get("structure_obj"),
                    which_bonds="all",
                )
            ],
            id=self.id("local-env-lobsterpy"),
        )

        return {
            "analysis-mode": analysis_mode,
            "local-envs": local_envs,
        }

    def layout(self):
        """Return the layout of the component."""
        # Get the sub-layouts
        # and create the main layout
        sub_layouts = self._sub_layouts

        controls = Columns(
            [
                Column(
                    [
                        sub_layouts["analysis-mode"],
                    ]
                )
            ]
        )

        # Create the local environments div
        local_envs_header = H4(
            "Local Environments identified via LobsterEnv",
            id=self.id("local-envs-text"),
            style={"display": "inline-block"},
        )
        local_envs_div = Columns([Column([sub_layouts["local-envs"]])])

        return Column(
            [
                controls,
                html.Br(),
                local_envs_header,
                local_envs_div,
            ]
        )

    @staticmethod
    def _get_all_inputs(
        data: dict | None,
    ) -> dict:
        data = data or {}

        charge_obj = data.get("charge_obj")
        icohplist_obj = data.get("icohplist_obj")
        structure_obj = data.get("structure_obj")

        if charge_obj and isinstance(charge_obj, dict):
            data["charge_obj"] = Charge.from_dict(charge_obj)

        if icohplist_obj and isinstance(icohplist_obj, dict):
            data["icohplist_obj"] = Icohplist.from_dict(icohplist_obj)

        if structure_obj and isinstance(structure_obj, dict):
            data["structure_obj"] = Structure.from_dict(structure_obj)

        return data

    @staticmethod
    def get_lobster_local_envs(
        charge_obj, icohplist_obj, structure_obj, which_bonds="all"
    ) -> str:
        """Get text description of local environments

        Args:
            input_dict: Dictionary containing the pymatgen objects.

        Returns:
            A string describing the local environments.
        """
        # Get the local environments using LobsterEnv

        spacegroup_analyzer = SpacegroupAnalyzer(structure=structure_obj)
        symm_struct = spacegroup_analyzer.get_symmetrized_structure()
        inequivalent_indices = [
            indices[0] for indices in symm_struct.equivalent_indices
        ]
        wyckoffs = symm_struct.wyckoff_symbols

        if which_bonds == "all":
            additional_condition = 0
            only_cation_environments = False
        else:
            additional_condition = 1
            only_cation_environments = True

        additional_condition = 0 if which_bonds == "all" else 1

        if icohplist_obj.are_cobis or icohplist_obj.are_coops:
            noise_cutoff = 0.001
        else:
            noise_cutoff = 0.1

        chem_env = LobsterNeighbors(
            filename_icohp=None,
            obj_icohp=icohplist_obj,
            structure=structure_obj,
            additional_condition=additional_condition,
            perc_strength_icohp=0.1,
            filename_charge=None,
            obj_charge=charge_obj,
            valences_from_charges=True,
            adapt_extremum_to_add_cond=True,
            are_cobis=icohplist_obj.are_cobis,
            are_coops=icohplist_obj.are_coops,
            noise_cutoff=noise_cutoff,
            which_charge="Mulliken",
        )

        lse = chem_env.get_light_structure_environment(
            only_cation_environments=only_cation_environments
        )

        save_dir_name = f"{chem_env.structure.composition.reduced_formula}"

        envs = []  # list of local environments
        summary_table = []
        for site_ix, env in enumerate(lse.coordination_environments):
            if site_ix in inequivalent_indices and env[0]["ce_symbol"]:
                # if env[0]["ce_symbol"]:
                data_list = []
                site_str = unicodeify_species(
                    chem_env.structure[site_ix].species_string
                )

                wyckoff_pos = wyckoffs[inequivalent_indices.index(site_ix)]

                try:
                    data_list.extend(
                        [
                            ["Specie", site_str],
                            ["Site Label", site_str + str(site_ix + 1)],
                            ["Wyckoff Positions", wyckoff_pos],
                            [
                                "Environment",
                                Description._coordination_environment_to_text(
                                    env[0]["ce_symbol"]
                                ).capitalize(),
                            ],
                            ["IUPAC Symbol", env[0]["ce_symbol"]],
                            ["CSM", float(round(env[0]["csm"], 5))],
                        ]
                    )

                except KeyError:
                    data_list.extend(
                        [
                            ["Specie", site_str],
                            ["Site Label", site_str + str(site_ix + 1)],
                            ["Wyckoff Positions", wyckoff_pos],
                            [
                                "Environment",
                                Description._coordination_environment_to_text(
                                    env[0]["ce_symbol"]
                                ).capitalize(),
                            ],
                            ["IUPAC Symbol", env[0]["ce_symbol"]],
                            ["CSM", "NA"],
                        ]
                    )

                local_env_data = chem_env.get_nn_info(chem_env.structure, site_ix)

                neighbour_sites = [i["site"] for i in local_env_data]
                central_site = chem_env.structure[site_ix]
                neighbour_weights = [
                    i["edge_properties"]["ICOHP"] for i in local_env_data
                ]
                charges = [charge_obj.mulliken[site_ix]]
                charges.extend(
                    [charge_obj.mulliken[i["site_index"]] for i in local_env_data]
                )

                # Create a molecule object for the local environment
                # and add the charges as a site property
                mol = Molecule.from_sites([central_site, *neighbour_sites])
                mol = mol.get_centered_molecule()

                # Add the charges as a site property (hover text)
                mol = mol.add_site_property("charge", charges)
                os.makedirs(name=save_dir_name, exist_ok=True)
                mol.to_file(
                    f"{save_dir_name}/{save_dir_name}_site_{site_ix}_{which_bonds}.xyz"
                )

                mg = MoleculeGraph.with_empty_graph(
                    molecule=mol,
                    name="bond_strength",
                    edge_weight_name="ICOHP",
                    edge_weight_units="eV",
                )
                for i in range(1, len(mol)):
                    # Add the bond strength as an edge weight (hover text)
                    mg.add_edge(0, i, weight=neighbour_weights[i - 1])

                view = html.Div(
                    [
                        StructureMoleculeComponent(
                            struct_or_mol=mg,
                            disable_callbacks=False,
                            id=f"{chem_env.structure.composition.reduced_formula}_site_{site_ix}",
                            scene_settings={
                                "enableZoom": True,
                                "defaultZoom": 1,
                            },
                            show_export_button=False,
                        )._sub_layouts["struct"]
                    ],
                    style={"width": "300px", "height": "300px"},
                )

                summary_table.append(data_list)

                data_list.append(["Interactive", view])

                envs.append(get_table(rows=data_list))

        envs_grouped = [envs[i : i + 2] for i in range(0, len(envs), 2)]
        analysis_contents = [
            Columns(
                [
                    Column(
                        html.Div(
                            e, style={"display": "flex", "justifyContent": "center"}
                        ),
                    )
                    for e in env_group
                ]
            )
            for env_group in envs_grouped
        ]

        dict_rows = [
            {item[0]: item[1] for item in row if item[0] != "Interactive"}
            for row in summary_table
        ]
        pd.DataFrame(dict_rows).to_csv(
            f"{save_dir_name}/{save_dir_name}_{which_bonds}_summary.csv"
        )

        return html.Div([html.Div(analysis_contents), html.Br()])

    def generate_callbacks(self, app, cache) -> None:
        """Register callback functions for this component."""

        @app.callback(
            Output(self.id("local-env-lobsterpy"), "children"),
            Input(self.id(), "data"),
            Input(self.get_kwarg_id("analysis-mode"), "value"),
        )
        def update_local_envs(data, label_select):
            """Update the local environments using LobsterEnv."""
            data = self._get_all_inputs(data)

            return self.get_lobster_local_envs(
                charge_obj=data.get("charge_obj"),
                icohplist_obj=data.get("icohplist_obj"),
                structure_obj=data.get("structure_obj"),
                which_bonds=(
                    label_select if isinstance(label_select, str) else label_select[0]
                ),
            )
