from __future__ import annotations

import json
import math
import re
import tempfile
import unittest
from pathlib import Path

from tiamat_py.io_json import load_json_project
from tiamat_py.io_ascii import parse_ascii_dat
from tiamat_py.io_dna import load_dna_project, parse_dna_records
from tiamat_py.model import Base, TiamatProject
from tiamat_py.io_pdb import convert_pdb_to_dna_files, load_pdb_sugar_center_project
from tiamat_py.model import basepair_state
from tiamat_py.viewer import (
    _apply_orientation,
    _apply_rotation_to_project,
    _append_free_strand,
    _allowed_sequence_bases,
    _base_edit_options,
    _base_terminal_roles,
    _basepair_dash,
    _blank_constraints_values,
    _blank_custom_constraints_values,
    _build_created_structure,
    _clipboard_render_scale,
    _clear_selection_on_empty_click,
    _constraint_profile_for_molecule_option,
    CUSTOM_CONSTRAINT_COLUMNS,
    Camera,
    _default_constraints_values,
    _apply_trackball_rotation,
    _drag_rotation_matrix_planar,
    _drag_translation_vector,
    _display_recent_file_path,
    _grid_axes_for_mode,
    _grid_fit_zoom,
    _grid_values,
    _helix_indices,
    _mark_viewport_for_fit,
    _menu_toggle_mark,
    _paste_structure_clipboard,
    _project_center,
    _projected_indices_in_box,
    _restore_history_state,
    _selection_clipboard_payload,
    _selected_strand_for_indices,
    _selected_edge_sets,
    _selection_targets_for_mode,
    _selection_status_text,
    _selection_center,
    _snap_create_helix_endpoint,
    _normalize_sequence_input,
    _project_constraints,
    _project_custom_constraints,
    _resolve_free_strand_creation,
    _snapshot_history_state,
    _strand_end_arrow_style,
    _strand_end_arrow_segment,
    _estimate_created_base_count,
    _estimate_free_strand_base_count,
    _validate_constraints_values,
    _rotation_matrix_z,
    FreeStrandEndpoint,
    ViewportState,
    load_project,
    project_points,
    render_svg,
    rotate_point,
)


SAMPLE_DAT = """\
4
1 AAAA 00000000 BBBB CCCC 0.0 0.0 0.0
2 BBBB AAAA 00000000 DDDD 1.0 0.0 0.0
3 CCCC 00000000 DDDD AAAA 0.0 1.0 0.0
4 DDDD CCCC 00000000 BBBB 1.0 1.0 0.0
"""

SAMPLE_DNA_PATH = Path("/Users/wyssuser/Desktop/Tiamat2_design/SSCCH/Ludlum_v0/14_33/ludlum_14_33.dna")
SAMPLE_DNAJSON_PATH = Path("/Users/wyssuser/Desktop/Tiamat2_design/SSCCH/Ludlum_v0/14_33/output.dnajson")
SAMPLE_DNA_PATH_VARIANT = Path("/Users/wyssuser/Desktop/Tiamat2_design/test.dna")
SAMPLE_DNAJSON_PATH_VARIANT = Path("/Users/wyssuser/Desktop/Tiamat2_design/test.dnajson")
SAMPLE_DNA_PATH_STANDALONE = Path("/Users/wyssuser/Desktop/Tiamat2_design/test.dna")
SAMPLE_DNAJSON_PATH_STANDALONE = Path("/Users/wyssuser/Desktop/Tiamat2_design/test.dnajson")
SAMPLE_RNA_DNA_PATH = Path("/Users/wyssuser/Desktop/Tiamat2_design/RNA_crisscross/V4_6565/V4_6565.dna")
SAMPLE_RNA_DNAJSON_PATH = Path("/Users/wyssuser/Desktop/Tiamat2_design/RNA_crisscross/V4_6565/V4_6565.dnajson")
SAMPLE_PDB_PATH = Path("/Users/wyssuser/Desktop/Tiamat2_design/ChemeraX_model/DNA.pdb")


class TiamatPyTests(unittest.TestCase):
    def test_parse_ascii_dat_builds_expected_graph(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        summary = project.summary()
        self.assertEqual(summary["bases"], 4)
        self.assertEqual(summary["strands"], 2)
        self.assertEqual(project.bases[1].down, 2)
        self.assertEqual(project.bases[1].across, 3)
        self.assertEqual(project.strands[0].base_indices, [1, 2])
        self.assertEqual(project.strands[1].base_indices, [3, 4])

    def test_viewer_loader_rejects_dat_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.dat"
            path.write_text(SAMPLE_DAT, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, r"\.dna, \.dnajson, or \.pdb"):
                load_project(path)

    def test_strand_order_uses_base_with_no_up_as_start(self) -> None:
        project = TiamatProject(
            bases={
                10: Base(index=10, object_id="10", x=1.0, y=0.0, z=0.0, up=20, down=30),
                20: Base(index=20, object_id="20", x=0.0, y=0.0, z=0.0, up=None, down=10),
                30: Base(index=30, object_id="30", x=2.0, y=0.0, z=0.0, up=10, down=None),
            }
        )
        self.assertEqual(project.strands[0].base_indices, [20, 10, 30])

    def test_circular_strand_uses_smallest_index_for_label_anchor(self) -> None:
        project = TiamatProject(
            bases={
                10: Base(index=10, object_id="10", x=0.0, y=0.0, z=0.0, up=30, down=20),
                20: Base(index=20, object_id="20", x=1.0, y=0.0, z=0.0, up=10, down=30),
                30: Base(index=30, object_id="30", x=2.0, y=0.0, z=0.0, up=20, down=10),
            }
        )
        self.assertEqual(project.strands[0].base_indices, [10, 20, 30])

    def test_delete_bases_breaks_neighbors_and_clears_across(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, up=None, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=3, across=4),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, up=2, down=None),
                4: Base(index=4, object_id="4", x=1.0, y=1.0, z=0.0, across=2),
            }
        )

        project.delete_bases({2})

        self.assertNotIn(2, project.bases)
        self.assertIsNone(project.bases[1].down)
        self.assertIsNone(project.bases[3].up)
        self.assertIsNone(project.bases[4].across)
        self.assertEqual([strand.base_indices for strand in project.strands], [[1], [3], [4]])

    def test_delete_bases_can_empty_project(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0),
            }
        )

        project.delete_bases({1})

        self.assertEqual(project.bases, {})
        self.assertEqual(project.strands, [])
        self.assertEqual(project.summary()["bases"], 0)

    def test_set_nucleotide_updates_canonical_pair_and_n_clears_both(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, molecule="DNA", across=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, molecule="DNA", across=1),
            }
        )
        project.set_nucleotide(1, "A")
        self.assertEqual(project.bases[1].nucleotide, "A")
        self.assertEqual(project.bases[2].nucleotide, "T")
        self.assertEqual(basepair_state(project.bases[1], project.bases[2]), "canonical")

        project.set_nucleotide(1, None)
        self.assertIsNone(project.bases[1].nucleotide)
        self.assertIsNone(project.bases[2].nucleotide)
        self.assertEqual(basepair_state(project.bases[1], project.bases[2]), "undefined")

    def test_set_nucleotide_preserves_rna_gu_wobble_pair(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, molecule="RNA", across=2, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, molecule="RNA", across=1, nucleotide="U"),
            }
        )
        project.set_nucleotide(1, "G")
        self.assertEqual(project.bases[1].nucleotide, "G")
        self.assertEqual(project.bases[2].nucleotide, "U")
        self.assertEqual(basepair_state(project.bases[1], project.bases[2]), "wobble")
        self.assertEqual(_basepair_dash(project, 1, 2), (4, 2))

    def test_create_and_delete_complementary_follow_selection_rules(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, molecule="RNA", nucleotide="G"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, molecule="RNA", nucleotide="U"),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, molecule="DNA", nucleotide="A"),
                4: Base(index=4, object_id="4", x=3.0, y=0.0, z=0.0, molecule="DNA", nucleotide="C"),
            }
        )
        self.assertTrue(project.can_create_complementary({1, 2}))
        self.assertTrue(project.create_complementary({1, 2}))
        self.assertEqual(project.bases[1].across, 2)
        self.assertEqual(project.bases[2].across, 1)
        self.assertTrue(project.can_delete_complementary({1, 2}))
        self.assertTrue(project.delete_complementary({1, 2}))
        self.assertIsNone(project.bases[1].across)
        self.assertIsNone(project.bases[2].across)
        self.assertFalse(project.can_create_complementary({3, 4}))

    def test_create_and_delete_down_update_up_down_links(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, up=None, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=None),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, up=None, down=None),
            }
        )
        self.assertTrue(project.can_create_down({2, 3}))
        self.assertTrue(project.create_down({2, 3}))
        self.assertEqual(project.bases[2].down, 3)
        self.assertEqual(project.bases[3].up, 2)
        self.assertTrue(project.can_delete_down({2}))
        self.assertTrue(project.delete_down({2}))
        self.assertIsNone(project.bases[2].down)
        self.assertIsNone(project.bases[3].up)

    def test_set_and_reset_backbone_color_only_affect_outgoing_backbone(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1),
            }
        )
        changed = project.set_backbone_color({1, 2}, "#123456")
        self.assertEqual(changed, 1)
        self.assertEqual(project.bases[1].backbone_color, "#123456")
        self.assertIsNone(project.bases[2].backbone_color)

        changed = project.reset_backbone_color({1, 2})
        self.assertEqual(changed, 1)
        self.assertIsNone(project.bases[1].backbone_color)

    def test_create_and_delete_sticky_end_follow_selection_rules(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=3, nucleotide="T"),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, up=2, down=4, nucleotide="C"),
                4: Base(index=4, object_id="4", x=3.0, y=0.0, z=0.0, up=3, nucleotide="G"),
                5: Base(index=5, object_id="5", x=0.0, y=1.0, z=0.0, down=6, nucleotide="C"),
                6: Base(index=6, object_id="6", x=1.0, y=1.0, z=0.0, up=5, down=7, nucleotide="G"),
                7: Base(index=7, object_id="7", x=2.0, y=1.0, z=0.0, up=6, down=8, nucleotide="A"),
                8: Base(index=8, object_id="8", x=3.0, y=1.0, z=0.0, up=7, nucleotide="T"),
            }
        )
        selection = {1, 2, 3, 4, 5, 6, 7, 8}
        self.assertTrue(project.can_create_sticky_end(selection))
        created, error = project.create_sticky_end(selection)
        self.assertTrue(created)
        self.assertIsNone(error)
        groups = project.sticky_end_groups()
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["left_indices"], [1, 2, 3, 4])
        self.assertEqual(groups[0]["right_indices"], [5, 6, 7, 8])
        self.assertTrue(project.can_delete_sticky_end({1}))
        self.assertTrue(project.delete_sticky_end({1}))
        self.assertEqual(project.sticky_end_groups(), [])

    def test_create_sticky_end_reports_non_complementary_segments(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, nucleotide="A"),
                3: Base(index=3, object_id="3", x=0.0, y=1.0, z=0.0, down=4, nucleotide="A"),
                4: Base(index=4, object_id="4", x=1.0, y=1.0, z=0.0, up=3, nucleotide="A"),
            }
        )
        created, error = project.create_sticky_end({1, 2, 3, 4})
        self.assertFalse(created)
        self.assertEqual(error, "Selected sequences are not reverse complementary.")
        self.assertEqual(project.sticky_end_groups(), [])

    def test_backbone_color_round_trips_through_history_snapshot(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1),
            }
        )
        project.set_backbone_color({1}, "#abcdef")

        snapshot = _snapshot_history_state(project, {1}, 1)
        restored_project, selected_indices, primary_selected = _restore_history_state(snapshot)

        self.assertIsNotNone(restored_project)
        self.assertEqual(restored_project.bases[1].backbone_color, "#abcdef")
        self.assertEqual(selected_indices, {1})
        self.assertEqual(primary_selected, 1)

    def test_fill_unpaired_generic_bases_only_changes_unpaired_n_and_maps_tu(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, molecule="DNA", nucleotide=None, across=None),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, molecule="RNA", nucleotide=None, across=None),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, molecule="DNA", nucleotide="A", across=None),
                4: Base(index=4, object_id="4", x=3.0, y=0.0, z=0.0, molecule="DNA", nucleotide=None, across=5),
                5: Base(index=5, object_id="5", x=4.0, y=0.0, z=0.0, molecule="DNA", nucleotide=None, across=4),
            }
        )
        changed = project.fill_unpaired_generic_bases("T")
        self.assertEqual(changed, 2)
        self.assertEqual(project.bases[1].nucleotide, "T")
        self.assertEqual(project.bases[2].nucleotide, "U")
        self.assertEqual(project.bases[3].nucleotide, "A")
        self.assertIsNone(project.bases[4].nucleotide)
        self.assertIsNone(project.bases[5].nucleotide)

    def test_selection_clipboard_payload_keeps_only_internal_links(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2, across=4, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=3, backbone_color="#abcdef", nucleotide="T"),
                3: Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, up=2, nucleotide="C"),
                4: Base(index=4, object_id="4", x=0.0, y=1.0, z=0.0, across=1, nucleotide="T"),
            }
        )
        payload = _selection_clipboard_payload(project, {1, 2})
        self.assertIsNotNone(payload)
        bases_payload = payload["bases"]
        self.assertEqual(len(bases_payload), 2)
        self.assertIsNone(bases_payload[0]["across"])
        self.assertEqual(bases_payload[0]["down"], 2)
        self.assertIsNone(bases_payload[1]["down"])
        self.assertIsNone(bases_payload[1]["backbone_color"])

    def test_paste_structure_clipboard_duplicates_selection_in_place(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, across=3, nucleotide="T"),
                3: Base(index=3, object_id="3", x=1.0, y=1.0, z=0.0, across=2, nucleotide="A"),
            }
        )
        payload = _selection_clipboard_payload(project, {1, 2, 3})
        pasted_project, new_indices = _paste_structure_clipboard(project, payload, (0.0, 0.0, 0.0))
        self.assertEqual(len(new_indices), 3)
        new_first = pasted_project.bases[min(new_indices)]
        self.assertAlmostEqual(new_first.x, 0.0)
        self.assertAlmostEqual(new_first.y, 0.0)
        self.assertAlmostEqual(new_first.z, 0.0)
        pasted_down = pasted_project.bases[new_first.down]
        self.assertEqual(pasted_down.up, new_first.index)
        pasted_pair = pasted_project.bases[pasted_down.across]
        self.assertEqual(pasted_pair.across, pasted_down.index)

    def test_selection_targets_cover_base_pair_strand_and_helix(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, up=None, down=2, across=3),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=None),
                3: Base(index=3, object_id="3", x=0.0, y=1.0, z=0.0, up=None, down=4, across=1),
                4: Base(index=4, object_id="4", x=1.0, y=1.0, z=0.0, up=3, down=None, across=5),
                5: Base(index=5, object_id="5", x=0.0, y=2.0, z=0.0, up=None, down=6, across=4),
                6: Base(index=6, object_id="6", x=1.0, y=2.0, z=0.0, up=5, down=None),
            }
        )

        self.assertEqual(_selection_targets_for_mode(project, 1, "base"), {1})
        self.assertEqual(_selection_targets_for_mode(project, 1, "base_pair"), {1, 3})
        self.assertEqual(_selection_targets_for_mode(project, 1, "strand"), {1, 2})
        self.assertEqual(_helix_indices(project, 1), {1, 2, 3, 4, 5, 6})
        self.assertEqual(_selection_targets_for_mode(project, 1, "helix"), {1, 2, 3, 4, 5, 6})

    def test_projected_indices_in_box_uses_rectangle_bounds(self) -> None:
        projected = {
            1: (10.0, 10.0, 0.0),
            2: (25.0, 20.0, 0.0),
            3: (80.0, 70.0, 0.0),
        }
        self.assertEqual(_projected_indices_in_box(projected, (5, 5), (30, 25)), {1, 2})

    def test_drag_translation_vector_respects_planar_mode(self) -> None:
        translation = _drag_translation_vector("xy", Camera(rot_x=0.0, rot_y=0.0, rot_z=0.0), dx=10.0, dy=20.0, scale=10.0)
        self.assertEqual(translation, (1.0, -2.0, 0.0))

    def test_drag_translation_vector_can_constrain_to_single_axis(self) -> None:
        translation = _drag_translation_vector(
            "xy",
            Camera(rot_x=0.0, rot_y=0.0, rot_z=0.0),
            dx=10.0,
            dy=20.0,
            scale=10.0,
            constrain_axis=True,
        )
        self.assertEqual(translation, (0.0, -2.0, 0.0))

    def test_base_edit_options_follow_molecule_type(self) -> None:
        self.assertEqual(
            _base_edit_options("DNA"),
            (
                ("Adenine", "A"),
                ("Thymine", "T"),
                ("Cytosine", "C"),
                ("Guanine", "G"),
                ("Generic", "N"),
            ),
        )
        self.assertEqual(
            _base_edit_options("RNA"),
            (
                ("Adenine", "A"),
                ("Uracil", "U"),
                ("Cytosine", "C"),
                ("Guanine", "G"),
                ("Generic", "N"),
            ),
        )

    def test_allowed_sequence_bases_follow_molecule_type(self) -> None:
        self.assertEqual(_allowed_sequence_bases("DNA"), ("A", "T", "C", "G", "N"))
        self.assertEqual(_allowed_sequence_bases("RNA"), ("A", "U", "C", "G", "N"))

    def test_normalize_sequence_input_ignores_spaces_and_reports_invalid_chars(self) -> None:
        normalized, invalid = _normalize_sequence_input("act gcn\nA$?", _allowed_sequence_bases("DNA"))
        self.assertEqual(normalized, "ACTGCNA")
        self.assertEqual(invalid, ["$", "?"])

    def test_selected_strand_for_indices_requires_one_strand(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        strand = _selected_strand_for_indices(project, {1, 2})
        self.assertIsNotNone(strand)
        self.assertEqual(strand.base_indices, [1, 2])
        self.assertIsNone(_selected_strand_for_indices(project, {1, 3}))

    def test_empty_click_clears_selection_only_in_direct_select_modes(self) -> None:
        for mode in ("helix", "strand", "base_pair", "base"):
            self.assertTrue(_clear_selection_on_empty_click(mode))
        for mode in ("box", "move", "rotate"):
            self.assertFalse(_clear_selection_on_empty_click(mode))

    def test_base_terminal_roles_identifies_five_and_three_prime_ends(self) -> None:
        five = Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, up=None, down=2)
        three = Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=None)
        isolated = Base(index=3, object_id="3", x=2.0, y=0.0, z=0.0, up=None, down=None)
        internal = Base(index=4, object_id="4", x=3.0, y=0.0, z=0.0, up=3, down=5)
        self.assertEqual(_base_terminal_roles(five), frozenset({"5"}))
        self.assertEqual(_base_terminal_roles(three), frozenset({"3"}))
        self.assertEqual(_base_terminal_roles(isolated), frozenset({"5", "3"}))
        self.assertEqual(_base_terminal_roles(internal), frozenset())

    def test_resolve_free_strand_creation_orients_from_three_to_five_prime(self) -> None:
        start = FreeStrandEndpoint(point=(0.0, 0.0, 0.0), base_index=10, terminal_roles=frozenset({"3"}), molecule="RNA")
        end = FreeStrandEndpoint(point=(3.0, 0.0, 0.0), base_index=11, terminal_roles=frozenset({"5"}), molecule="RNA")
        plan = _resolve_free_strand_creation(start, end)
        self.assertEqual(plan["start_anchor"], 10)
        self.assertEqual(plan["end_anchor"], 11)
        self.assertEqual(plan["molecule"], "RNA")
        self.assertEqual(_estimate_free_strand_base_count(start.point, end.point), 5)

    def test_append_free_strand_links_new_segment_between_existing_ends(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, up=None, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, down=None),
                3: Base(index=3, object_id="3", x=5.0, y=0.0, z=0.0, up=None, down=4),
                4: Base(index=4, object_id="4", x=6.0, y=0.0, z=0.0, up=3, down=None),
            }
        )
        created = _append_free_strand(
            project,
            start_point=project.bases[2].position,
            end_point=project.bases[3].position,
            count=2,
            molecule="DNA",
            start_anchor=2,
            end_anchor=3,
        )
        self.assertEqual(created, [5, 6])
        self.assertEqual(project.bases[2].down, 5)
        self.assertEqual(project.bases[5].up, 2)
        self.assertEqual(project.bases[5].down, 6)
        self.assertEqual(project.bases[6].down, 3)
        self.assertEqual(project.bases[3].up, 6)
        self.assertEqual(len(project.strands), 1)

    def test_planar_rotation_uses_view_axes_for_xz_drag(self) -> None:
        horizontal_rotation = _drag_rotation_matrix_planar("xz", dx=90.0, dy=0.0, sensitivity=1.0)
        rotated = _apply_orientation(horizontal_rotation, 1.0, 0.0, 0.0)
        self.assertAlmostEqual(rotated[0], 0.0, places=6)
        self.assertAlmostEqual(rotated[1], 1.0, places=6)
        self.assertAlmostEqual(rotated[2], 0.0, places=6)

        vertical_rotation = _drag_rotation_matrix_planar("xz", dx=0.0, dy=-90.0, sensitivity=1.0)
        rotated = _apply_orientation(vertical_rotation, 0.0, 1.0, 0.0)
        self.assertAlmostEqual(rotated[0], 0.0, places=6)
        self.assertAlmostEqual(rotated[1], 0.0, places=6)
        self.assertAlmostEqual(rotated[2], 1.0, places=6)

    def test_planar_rotation_can_constrain_to_single_axis(self) -> None:
        constrained = _drag_rotation_matrix_planar("xz", dx=90.0, dy=-45.0, sensitivity=1.0, constrain_axis=True)
        rotated = _apply_orientation(constrained, 1.0, 0.0, 0.0)
        self.assertAlmostEqual(rotated[0], 0.0, places=6)
        self.assertAlmostEqual(rotated[1], 1.0, places=6)
        self.assertAlmostEqual(rotated[2], 0.0, places=6)

    def test_apply_rotation_to_project_rotates_about_selection_center(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=1.0, y=0.0, z=0.0),
                2: Base(index=2, object_id="2", x=-1.0, y=0.0, z=0.0),
            }
        )
        selected = {1, 2}
        center = _selection_center(project, selected)
        original_positions = {index: project.bases[index].position for index in selected}
        _apply_rotation_to_project(project, original_positions, center, _rotation_matrix_z(90.0))
        self.assertAlmostEqual(project.bases[1].x, 0.0, places=6)
        self.assertAlmostEqual(project.bases[1].y, 1.0, places=6)
        self.assertAlmostEqual(project.bases[2].x, 0.0, places=6)
        self.assertAlmostEqual(project.bases[2].y, -1.0, places=6)

    def test_sequence_generation_assigns_complements(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.generate_sequences(seed=1)
        self.assertIsNotNone(project.bases[1].nucleotide)
        self.assertEqual(project.bases[3].nucleotide, {"A": "T", "T": "A", "C": "G", "G": "C"}[project.bases[1].nucleotide])
        self.assertEqual(project.bases[4].nucleotide, {"A": "T", "T": "A", "C": "G", "G": "C"}[project.bases[2].nucleotide])

    def test_sequence_generation_respects_unique_window_and_repetition_limits(self) -> None:
        bases = {
            index: Base(
                index=index,
                object_id=str(index),
                x=float(index),
                y=0.0,
                z=0.0,
                molecule="DNA",
                up=index - 1 if index > 1 else None,
                down=index + 1 if index < 12 else None,
            )
            for index in range(1, 13)
        }
        project = TiamatProject(bases=bases)
        project.generate_sequences(
            seed=1,
            unique_sequence_limit=4,
            repetition_limit=3,
            gc_percentage=0.5,
            gc_tolerance=0.5,
            timeout=1.0,
        )
        sequence = project.strands[0].sequence(project)
        self.assertEqual(len(sequence), 12)
        for start in range(0, len(sequence) - 4 + 1):
            self.assertEqual(set(sequence[start : start + 4]), set("ATCG"))
        for start in range(0, len(sequence) - 3 + 1):
            self.assertNotEqual(len(set(sequence[start : start + 3])), 1)

    def test_sequence_generation_respects_gc_range(self) -> None:
        bases = {
            index: Base(
                index=index,
                object_id=str(index),
                x=float(index),
                y=0.0,
                z=0.0,
                molecule="DNA",
                up=index - 1 if index > 1 else None,
                down=index + 1 if index < 40 else None,
            )
            for index in range(1, 41)
        }
        project = TiamatProject(bases=bases)
        project.generate_sequences(
            seed=2,
            unique_sequence_limit=10,
            repetition_limit=6,
            gc_percentage=0.5,
            gc_tolerance=0.05,
            timeout=1.0,
        )
        sequence = project.strands[0].sequence(project)
        gc_ratio = sum(1 for base in sequence if base in {"G", "C"}) / len(sequence)
        self.assertGreaterEqual(gc_ratio, 0.45)
        self.assertLessEqual(gc_ratio, 0.55)

    def test_json_round_trip(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.generate_sequences(seed=2)
        loaded = TiamatProject.from_dict(json.loads(project.to_json()))
        self.assertEqual(loaded.export_sequences_text(), project.export_sequences_text())

    def test_default_and_blank_constraints_values_have_expected_shape(self) -> None:
        defaults = _default_constraints_values()
        blanks = _blank_constraints_values()
        custom_blanks = _blank_custom_constraints_values()
        self.assertEqual(defaults["rotation_per_base_pair_deg"]["b_median"], "-34.2857")
        self.assertEqual(defaults["helicity_bp_per_turn"]["b_median"], "10.50")
        self.assertEqual(defaults["diameter_nm"]["a_median"], "2.3")
        self.assertEqual(blanks["rotation_per_base_pair_deg"]["b_median"], "")
        self.assertEqual(blanks["helicity_bp_per_turn"]["b_median"], "")
        self.assertEqual(blanks["minor_groove_angle_deg"]["a_median"], "")
        self.assertEqual(custom_blanks["rotation_per_base_pair_deg"]["x_median"], "")
        self.assertEqual(custom_blanks["helicity_bp_per_turn"]["x_median"], "")
        self.assertEqual(custom_blanks["minor_groove_angle_deg"]["x_median"], "")
        self.assertNotIn("b_error", defaults["rotation_per_base_pair_deg"])
        self.assertNotIn("x_error", custom_blanks["rotation_per_base_pair_deg"])

    def test_validate_constraints_values_accepts_numeric_and_blank(self) -> None:
        values = _blank_constraints_values()
        values["rotation_per_base_pair_deg"]["b_median"] = "-34.2857"
        values["rise_per_base_pair_nm"]["a_median"] = "0.29"
        normalized = _validate_constraints_values(values)
        self.assertEqual(normalized["rotation_per_base_pair_deg"]["b_median"], "-34.2857")
        self.assertEqual(normalized["helicity_bp_per_turn"]["b_median"], "10.50")
        self.assertEqual(normalized["inclination_deg"]["a_median"], "")
        values["diameter_nm"]["b_median"] = "oops"
        with self.assertRaisesRegex(ValueError, "must be numeric or blank"):
            _validate_constraints_values(values)

        custom_values = _blank_custom_constraints_values()
        custom_values["helicity_bp_per_turn"]["x_median"] = "9.50"
        custom_normalized = _validate_constraints_values(custom_values, columns=CUSTOM_CONSTRAINT_COLUMNS)
        self.assertEqual(custom_normalized["helicity_bp_per_turn"]["x_median"], "9.50")
        self.assertEqual(custom_normalized["rotation_per_base_pair_deg"]["x_median"], "-37.8947")
        custom_values["helicity_bp_per_turn"]["x_median"] = "0"
        with self.assertRaisesRegex(ValueError, "non-zero"):
            _validate_constraints_values(custom_values, columns=CUSTOM_CONSTRAINT_COLUMNS)

    def test_project_constraints_loads_metadata_when_present(self) -> None:
        project = TiamatProject(
            bases={1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0)},
            metadata={
                "constraints": {
                    "diameter_nm": {
                        "b_median": "1.9",
                        "a_median": "2.2",
                    }
                }
            },
        )
        constraints = _project_constraints(project)
        self.assertEqual(constraints["diameter_nm"]["b_median"], "1.9")
        self.assertEqual(constraints["diameter_nm"]["a_median"], "2.2")
        self.assertEqual(constraints["helicity_bp_per_turn"]["b_median"], "10.50")
        self.assertNotIn("a_error", constraints["diameter_nm"])
        self.assertEqual(constraints["rotation_per_base_pair_deg"]["b_median"], "-34.2857")

    def test_project_custom_constraints_loads_metadata_when_present(self) -> None:
        project = TiamatProject(
            bases={1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0)},
            metadata={
                "custom_constraints": {
                    "diameter_nm": {
                        "x_median": "2.5",
                    }
                }
            },
        )
        constraints = _project_custom_constraints(project)
        self.assertEqual(constraints["diameter_nm"]["x_median"], "2.5")
        self.assertEqual(constraints["rotation_per_base_pair_deg"]["x_median"], "")
        self.assertEqual(constraints["helicity_bp_per_turn"]["x_median"], "")

    def test_estimate_created_base_count_uses_b_type_rise(self) -> None:
        self.assertEqual(_estimate_created_base_count(3.32, _default_constraints_values()), 10)
        self.assertEqual(_estimate_created_base_count(0.0, _default_constraints_values()), 60)

    def test_constraint_profile_for_molecule_option_maps_family_and_molecules(self) -> None:
        profile = _constraint_profile_for_molecule_option("DNA-RNA A", _default_constraints_values())
        self.assertEqual(profile["strand_molecules"], ("DNA", "RNA"))
        self.assertEqual(profile["helix_family"], "a")
        self.assertAlmostEqual(profile["rise"], 0.29)
        self.assertAlmostEqual(profile["inclination"], 19.0)
        self.assertAlmostEqual(profile["minor_groove_angle"], 106.17)
        self.assertNotIn("chord_length", profile)

    def test_constraint_profile_for_molecule_option_supports_x_type(self) -> None:
        custom = _blank_custom_constraints_values()
        custom["rotation_per_base_pair_deg"]["x_median"] = "-40"
        custom["rise_per_base_pair_nm"]["x_median"] = "0.4"
        custom["diameter_nm"]["x_median"] = "2.1"
        custom["inclination_deg"]["x_median"] = "12"
        custom["minor_groove_angle_deg"]["x_median"] = "120"
        profile = _constraint_profile_for_molecule_option("RNA-DNA X", custom)
        self.assertEqual(profile["strand_molecules"], ("RNA", "DNA"))
        self.assertEqual(profile["helix_family"], "x")
        self.assertAlmostEqual(profile["rise"], 0.4)
        self.assertAlmostEqual(profile["minor_groove_angle"], 120.0)

    def test_build_created_structure_creates_antiparallel_helix(self) -> None:
        created = _build_created_structure(
            start_point=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            count=3,
            backbone_rotation=0.0,
            molecule_option="DNA-DNA B",
            structure_type="Helix",
            constraints_values=_default_constraints_values(),
            start_index=10,
        )
        self.assertEqual(sorted(created), [10, 11, 12, 13, 14, 15])
        self.assertEqual(created[10].molecule, "DNA")
        self.assertEqual(created[11].molecule, "DNA")
        self.assertEqual(created[10].down, 12)
        self.assertIsNone(created[10].up)
        self.assertEqual(created[12].up, 10)
        self.assertEqual(created[11].down, None)
        self.assertEqual(created[11].up, 13)
        self.assertEqual(created[10].across, 11)
        self.assertEqual(created[12].across, 13)
        self.assertEqual(created[14].across, 15)
        expected_pair_distance = 2.0 * math.sin(math.radians(135.92 / 2.0))
        self.assertAlmostEqual(math.dist(created[10].position, created[11].position), expected_pair_distance, places=3)
        first_angle = math.atan2(created[10].z, created[10].y)
        second_angle = math.atan2(created[12].z, created[12].y)
        delta = ((second_angle - first_angle + math.pi) % (2.0 * math.pi)) - math.pi
        self.assertGreater(delta, 0.0)
        pair_angle = math.atan2(created[11].z, created[11].y)
        pair_delta = ((pair_angle - first_angle + math.pi) % (2.0 * math.pi)) - math.pi
        self.assertAlmostEqual(math.degrees(pair_delta), 135.92, delta=0.2)
        self.assertIsNotNone(created[10].plane_u)
        self.assertIsNotNone(created[10].plane_v)

    def test_build_created_structure_creates_single_strand_from_first_molecule(self) -> None:
        created = _build_created_structure(
            start_point=(0.0, 0.0, 0.0),
            direction=(0.0, 0.0, 0.0),
            count=2,
            backbone_rotation=30.0,
            molecule_option="RNA-DNA A",
            structure_type="Strand",
            constraints_values=_default_constraints_values(),
            start_index=1,
        )
        self.assertEqual(sorted(created), [1, 2])
        self.assertEqual(created[1].molecule, "RNA")
        self.assertEqual(created[2].molecule, "RNA")
        self.assertEqual(created[1].down, 2)
        self.assertIsNone(created[1].up)
        self.assertIsNone(created[1].across)
        self.assertIsNotNone(created[1].plane_u)
        self.assertIsNotNone(created[1].plane_v)

    def test_build_created_structure_applies_a_form_inclination_to_base_plane(self) -> None:
        created = _build_created_structure(
            start_point=(0.0, 0.0, 0.0),
            direction=(1.0, 0.0, 0.0),
            count=2,
            backbone_rotation=0.0,
            molecule_option="RNA-RNA A",
            structure_type="Helix",
            constraints_values=_default_constraints_values(),
            start_index=20,
        )
        base = created[20]
        paired_base = created[21]
        self.assertIsNotNone(base.plane_u)
        self.assertIsNotNone(base.plane_v)
        expected_offset = 1.15 * math.tan(math.radians(19.0))
        self.assertAlmostEqual(base.x, expected_offset, places=3)
        self.assertAlmostEqual(paired_base.x, -expected_offset, places=3)
        normal = (
            base.plane_u[1] * base.plane_v[2] - base.plane_u[2] * base.plane_v[1],
            base.plane_u[2] * base.plane_v[0] - base.plane_u[0] * base.plane_v[2],
            base.plane_u[0] * base.plane_v[1] - base.plane_u[1] * base.plane_v[0],
        )
        normal_length = math.sqrt(sum(value * value for value in normal))
        unit_normal = tuple(value / normal_length for value in normal)
        angle = math.degrees(math.acos(max(-1.0, min(1.0, unit_normal[0]))))
        self.assertAlmostEqual(angle, 19.0, delta=0.6)

    def test_history_snapshot_round_trip_restores_project_and_selection(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.bases[1].nucleotide = "A"
        snapshot = _snapshot_history_state(project, {1, 3}, 3)
        restored_project, restored_selected, restored_primary = _restore_history_state(snapshot)
        self.assertIsNotNone(restored_project)
        self.assertEqual(restored_project.bases[1].nucleotide, "A")
        self.assertEqual(restored_project.summary(), project.summary())
        self.assertEqual(restored_selected, {1, 3})
        self.assertEqual(restored_primary, 3)

    def test_svg_render_contains_svg_root(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        svg = render_svg(project, Camera())
        self.assertIn("<svg", svg)
        self.assertIn("<circle", svg)
        self.assertIn("<polygon", svg)
        self.assertIn('<rect x="0" y="0" width="1200" height="900"', svg)

    def test_svg_render_supports_planar_modes(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        for mode in ("xy", "yz", "xz"):
            svg = render_svg(project, Camera(), mode=mode)
            self.assertIn("<svg", svg)
            self.assertIn("<circle", svg)

    def test_projection_supports_all_modes(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        for mode in ("xy", "yz", "xz", "3d"):
            projected = project_points(project, Camera(), 400, 300, mode=mode)
            self.assertEqual(set(projected), {1, 2, 3, 4})

    def test_snap_create_helix_endpoint_locks_to_nearest_axis(self) -> None:
        snapped = _snap_create_helix_endpoint(
            start_canvas=(200, 150),
            end_canvas=(260, 120),
            camera=Camera(),
            width=400,
            height=300,
            mode="xy",
            center=(0.0, 0.0, 0.0),
            scale=10.0,
        )
        self.assertEqual(snapped[1], 150)

    def test_projection_center_can_stay_fixed_across_geometry_changes(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0),
                2: Base(index=2, object_id="2", x=10.0, y=0.0, z=0.0),
            }
        )
        fixed_center = _project_center(project)
        before = project_points(project, Camera(zoom=10.0), 400, 300, mode="xy", center=fixed_center)
        project.bases[2].x = 20.0
        after = project_points(project, Camera(zoom=10.0), 400, 300, mode="xy", center=fixed_center)
        self.assertAlmostEqual(before[1][0], after[1][0], places=6)
        self.assertAlmostEqual(before[1][1], after[1][1], places=6)

    def test_mark_viewport_for_fit_can_reset_stale_pan(self) -> None:
        state = ViewportState(
            key="xy",
            title="XY",
            mode="xy",
            canvas=None,
            camera=Camera(zoom=10.0, pan_x=185.0, pan_y=-92.0),
            needs_fit=False,
            fit_zoom=10.0,
            scene_center=(3.0, 4.0, 5.0),
        )

        _mark_viewport_for_fit(state, recenter=True)

        self.assertTrue(state.needs_fit)
        self.assertEqual(state.fit_zoom, 0.0)
        self.assertIsNone(state.scene_center)
        self.assertEqual(state.camera.pan_x, 0.0)
        self.assertEqual(state.camera.pan_y, 0.0)

    def test_strand_end_arrow_segment_extends_from_three_prime_end(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        projected = project_points(project, Camera(), 400, 300, mode="xy")
        segment = _strand_end_arrow_segment(project, project.strands[0], projected)
        self.assertIsNotNone(segment)
        x1, _y1, x2, _y2 = segment
        self.assertGreater(x2, x1)

    def test_strand_end_arrow_style_scales_with_radius(self) -> None:
        small_extension, small_arrowshape, small_svg_length, small_svg_half_width = _strand_end_arrow_style(2.0)
        large_extension, large_arrowshape, large_svg_length, large_svg_half_width = _strand_end_arrow_style(8.0)
        self.assertGreater(large_extension, small_extension)
        self.assertGreater(large_arrowshape[0], small_arrowshape[0])
        self.assertGreater(large_svg_length, small_svg_length)
        self.assertGreater(large_svg_half_width, small_svg_half_width)

    def test_trackball_rotation_moves_x_axis_out_of_flat_plane(self) -> None:
        camera = Camera()
        rotated_camera = Camera(orientation=_apply_trackball_rotation(camera, dx=0.0, dy=40.0))
        x_axis = rotate_point(1.0, 0.0, 0.0, rotated_camera)
        self.assertNotAlmostEqual(x_axis[1], 0.0, places=6)

    def test_external_dnajson_schema_loads(self) -> None:
        payload = {
            "bases": [
                {
                    "id": 0,
                    "position": [0.0, 0.0, 0.0],
                    "molecule": "DNA",
                    "type": "Adenine",
                    "across": 1,
                    "up": None,
                    "down": None,
                },
                {
                    "id": 1,
                    "position": [1.0, 0.0, 0.0],
                    "molecule": "RNA",
                    "type": "Uracil",
                    "across": 0,
                    "up": None,
                    "down": None,
                },
            ]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.dnajson"
            path.write_text(json.dumps(payload), encoding="utf-8")
            project = load_json_project(path)
        self.assertEqual(project.bases[0].nucleotide, "A")
        self.assertEqual(project.bases[1].nucleotide, "U")
        self.assertEqual(project.bases[1].molecule, "RNA")

    def test_pdb_sugar_center_loader_builds_expected_duplex(self) -> None:
        project = load_pdb_sugar_center_project(SAMPLE_PDB_PATH)
        summary = project.summary()
        self.assertEqual(summary["bases"], 92)
        self.assertEqual(summary["strands"], 2)
        self.assertEqual(summary["paired_bases"], 92)
        self.assertEqual(project.strands[0].sequence(project), "ATCGATCGAGTCTAGCTGCATCTGCTAGCGACTAGCTGACTGTAGT")
        self.assertAlmostEqual(project.bases[1].x, 0.57986, places=4)
        self.assertAlmostEqual(project.bases[1].y, 6.81130, places=4)
        self.assertAlmostEqual(project.bases[1].z, 1.46572, places=4)

    def test_viewer_loader_accepts_pdb_files(self) -> None:
        project = load_project(SAMPLE_PDB_PATH)
        self.assertEqual(project.summary()["bases"], 92)
        self.assertEqual(project.bases[1].molecule, "DNA")

    def test_pdb_converter_writes_minimal_dna_that_round_trips(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            dna_path = Path(tmpdir) / "dna_from_pdb.dna"
            _project, output_path, sidecar_path = convert_pdb_to_dna_files(SAMPLE_PDB_PATH, dna_path=dna_path)
            self.assertEqual(output_path, dna_path)
            self.assertIsNotNone(sidecar_path)
            loaded = load_dna_project(output_path)
            summary = loaded.summary()
            self.assertEqual(summary["bases"], 92)
            self.assertEqual(summary["strands"], 2)
            self.assertEqual(summary["paired_bases"], 92)

    def test_svg_uses_rect_for_rna(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.bases[1].molecule = "RNA"
        svg = render_svg(project, Camera(), selected_indices={1})
        self.assertIn("<rect", svg)
        self.assertIn("<polygon", svg)

    def test_svg_uses_dashed_line_for_rna_wobble_pair(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, molecule="RNA", across=2, nucleotide="G"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, molecule="RNA", across=1, nucleotide="U"),
            }
        )
        svg = render_svg(project, Camera())
        self.assertIn('stroke-dasharray="4,2"', svg)

    def test_selected_edges_only_include_next_backbone_and_pair(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        selected_along, selected_across = _selected_edge_sets(project, {1, 2})
        self.assertEqual(selected_along, [(1, 2)])
        self.assertEqual(selected_across, [(1, 3), (2, 4)])

    def test_selection_status_shows_single_id(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        self.assertEqual(_selection_status_text(project, {2}, 2), "Selected base ID: 2")

    def test_selection_status_shows_count_and_distance_for_two(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        self.assertEqual(_selection_status_text(project, {1, 2}, 2), "Selected bases: 2 | Distance: 1.00 nm")

    def test_selection_status_shows_count_for_many(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        self.assertEqual(_selection_status_text(project, {1, 2, 3}, 3), "Selected bases: 3")

    def test_menu_toggle_mark_shows_check_only_when_enabled(self) -> None:
        self.assertEqual(_menu_toggle_mark(True), "\u2713")
        self.assertEqual(_menu_toggle_mark(False), "")

    def test_display_recent_file_path_compacts_long_paths(self) -> None:
        compact = _display_recent_file_path("/very/long/path/to/a/project/example_structure_file.dna", max_length=24)
        self.assertIn("...", compact)
        self.assertLessEqual(len(compact), 24)

    def test_grid_axes_follow_planar_view(self) -> None:
        self.assertEqual(_grid_axes_for_mode("xy"), ("x", "y"))
        self.assertEqual(_grid_axes_for_mode("yz"), ("y", "z"))
        self.assertEqual(_grid_axes_for_mode("xz"), ("x", "z"))

    def test_grid_values_include_zero_and_spacing_steps(self) -> None:
        self.assertEqual(_grid_values(-2.0, 2.0, 1.0), [-2.0, -1.0, 0.0, 1.0, 2.0])

    def test_grid_fit_zoom_uses_axis_span(self) -> None:
        zoom_xy = _grid_fit_zoom(
            {"x_min": -5.0, "x_max": 5.0, "y_min": -2.0, "y_max": 2.0, "z_min": -1.0, "z_max": 1.0, "spacing": 1.0},
            1000,
            600,
            "xy",
        )
        zoom_yz = _grid_fit_zoom(
            {"x_min": -5.0, "x_max": 5.0, "y_min": -2.0, "y_max": 2.0, "z_min": -1.0, "z_max": 1.0, "spacing": 1.0},
            1000,
            600,
            "yz",
        )
        self.assertGreater(zoom_yz, zoom_xy)

    def test_clipboard_render_scale_upsizes_small_views_and_keeps_large_views(self) -> None:
        self.assertEqual(_clipboard_render_scale(1200, 900, max_dimension=2400.0), 2.0)
        self.assertEqual(_clipboard_render_scale(3200, 1800, max_dimension=2400.0), 1.0)

    def test_svg_selection_keeps_backbone_color_and_basepair_width(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        svg = render_svg(project, Camera(), selected_indices={1})
        self.assertIn('stroke="#ffffff" stroke-width="1"', svg)
        self.assertIn('stroke="#ffffff" stroke-width="4"', svg)
        self.assertIn(project.strands[0].color, svg)

    def test_svg_honors_custom_display_sizes(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        svg = render_svg(project, Camera(), sphere_size=10.0, backbone_width=5, basepair_width=4)
        self.assertIn('stroke-width="5"', svg)
        self.assertIn('stroke-width="4"', svg)
        self.assertIn('r="10.00"', svg)

    def test_svg_bright_mode_uses_white_background_and_black_selection_outline(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        svg = render_svg(project, Camera(), selected_indices={1}, color_scheme="bright")
        self.assertIn('fill="#ffffff"', svg)
        self.assertIn('stroke="#111111" stroke-width="1"', svg)
        self.assertIn('stroke="#111111" stroke-width="4"', svg)

    def test_svg_uses_custom_backbone_color_for_segment_and_three_prime_arrow(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1),
            }
        )
        default_color = project.strands[0].color
        project.set_backbone_color({1}, "#123456")

        svg = render_svg(project, Camera())

        self.assertIn('stroke="#123456"', svg)
        self.assertNotIn(default_color, svg)

    def test_svg_show_sticky_end_toggle_controls_dashed_lines_and_boxes(self) -> None:
        project = TiamatProject(
            bases={
                1: Base(index=1, object_id="1", x=0.0, y=0.0, z=0.0, down=2, nucleotide="A"),
                2: Base(index=2, object_id="2", x=1.0, y=0.0, z=0.0, up=1, nucleotide="T"),
                3: Base(index=3, object_id="3", x=0.0, y=1.0, z=0.0, down=4, nucleotide="A"),
                4: Base(index=4, object_id="4", x=1.0, y=1.0, z=0.0, up=3, nucleotide="T"),
            }
        )
        created, error = project.create_sticky_end({1, 2, 3, 4})
        self.assertTrue(created)
        self.assertIsNone(error)
        sticky_color = project.sticky_end_groups()[0]["color"]

        shown_svg = render_svg(project, Camera(), show_sticky_end=True)
        hidden_svg = render_svg(project, Camera(), show_sticky_end=False)

        self.assertIn('stroke-dasharray="12,8"', shown_svg)
        self.assertIn(sticky_color, shown_svg)
        self.assertNotIn('stroke-dasharray="12,8"', hidden_svg)
        self.assertNotIn(sticky_color, hidden_svg)

    def test_rna_marker_scales_with_custom_size(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.bases[1].molecule = "RNA"
        svg = render_svg(project, Camera(), sphere_size=10.0)
        self.assertIn('width="20.00"', svg)

    def test_svg_scales_dna_marker_size_with_zoom(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        close_svg = render_svg(project, Camera(zoom=25.0), sphere_size=10.0)
        detail_svg = render_svg(project, Camera(zoom=250.0), sphere_size=10.0)
        close_radius = float(re.search(r'r="([0-9.]+)"', close_svg).group(1))
        detail_radius = float(re.search(r'r="([0-9.]+)"', detail_svg).group(1))
        self.assertLess(close_radius, detail_radius)

    def test_svg_scales_rna_marker_size_with_zoom(self) -> None:
        project = parse_ascii_dat(SAMPLE_DAT)
        project.bases[1].molecule = "RNA"
        close_svg = render_svg(project, Camera(zoom=25.0), sphere_size=10.0)
        detail_svg = render_svg(project, Camera(zoom=250.0), sphere_size=10.0)
        close_width = float(
            re.findall(
                r'<rect x="[^"]+" y="[^"]+" width="([0-9.]+)" height="[0-9.]+" fill="#[0-9a-fA-F]{6}" stroke="#111827"',
                close_svg,
            )[0]
        )
        detail_width = float(
            re.findall(
                r'<rect x="[^"]+" y="[^"]+" width="([0-9.]+)" height="[0-9.]+" fill="#[0-9a-fA-F]{6}" stroke="#111827"',
                detail_svg,
            )[0]
        )
        self.assertLess(close_width, detail_width)

    def test_dna_binary_records_parse_from_sample(self) -> None:
        if not SAMPLE_DNA_PATH.exists():
            self.skipTest("sample .dna fixture is not available")
        records = parse_dna_records(SAMPLE_DNA_PATH)
        self.assertEqual(len(records), 159)
        self.assertEqual(records[0].nucleotide, "C")
        self.assertEqual(records[-1].nucleotide, "T")

    def test_dna_binary_load_matches_sidecar_summary(self) -> None:
        if not SAMPLE_DNA_PATH.exists() or not SAMPLE_DNAJSON_PATH.exists():
            self.skipTest("sample .dna/.dnajson fixtures are not available")
        project = load_dna_project(SAMPLE_DNA_PATH)
        sidecar = load_json_project(SAMPLE_DNAJSON_PATH)
        self.assertEqual(project.summary(), sidecar.summary())
        self.assertEqual(project.metadata.get("source_format"), "dna")
        self.assertEqual(project.metadata.get("dna_sidecar_path"), str(SAMPLE_DNAJSON_PATH))

    def test_dna_binary_loads_generic_rna_variant(self) -> None:
        if not SAMPLE_DNA_PATH_VARIANT.exists() or not SAMPLE_DNAJSON_PATH_VARIANT.exists():
            self.skipTest("variant .dna/.dnajson fixtures are not available")
        records = parse_dna_records(SAMPLE_DNA_PATH_VARIANT)
        project = load_dna_project(SAMPLE_DNA_PATH_VARIANT)
        sidecar = load_json_project(SAMPLE_DNAJSON_PATH_VARIANT)
        self.assertGreater(len(records), 0)
        self.assertIn(records[0].nucleotide, {"A", "C", "G", "T", "U", None})
        self.assertIn(records[0].molecule, {"DNA", "RNA", None})
        self.assertTrue(any(record.molecule == "RNA" or record.nucleotide == "U" for record in records))
        self.assertEqual(project.summary(), sidecar.summary())
        self.assertEqual(project.metadata.get("dna_record_count"), len(records))
        self.assertEqual(project.metadata.get("source_format"), "dna")
        self.assertEqual(project.metadata.get("dna_sidecar_path"), str(SAMPLE_DNAJSON_PATH_VARIANT))

    def test_dna_binary_standalone_fallback_matches_record_order_topology(self) -> None:
        if not SAMPLE_DNA_PATH_STANDALONE.exists() or not SAMPLE_DNAJSON_PATH_STANDALONE.exists():
            self.skipTest("standalone .dna/.dnajson fixtures are not available")

        sidecar = load_json_project(SAMPLE_DNAJSON_PATH_STANDALONE)
        with tempfile.TemporaryDirectory() as tmpdir:
            isolated_dna = Path(tmpdir) / SAMPLE_DNA_PATH_STANDALONE.name
            isolated_dna.write_bytes(SAMPLE_DNA_PATH_STANDALONE.read_bytes())
            project = load_dna_project(isolated_dna)

        self.assertEqual(project.summary(), sidecar.summary())
        self.assertEqual(project.metadata.get("dna_topology"), "record_order_fallback")
        self.assertIsNone(project.metadata.get("dna_sidecar_path"))
        coordinate_map = {}
        unmatched = dict(sidecar.bases)
        for index, base in project.bases.items():
            best_index = None
            best_distance = None
            for candidate_index, candidate in unmatched.items():
                distance = math.dist(base.position, candidate.position)
                if distance > 1e-2:
                    continue
                if best_distance is None or distance < best_distance:
                    best_distance = distance
                    best_index = candidate_index
            self.assertIsNotNone(best_index)
            coordinate_map[index] = best_index
            unmatched.pop(best_index)

        for index, sidecar_index in coordinate_map.items():
            loaded_base = project.bases[index]
            sidecar_base = sidecar.bases[sidecar_index]
            self.assertEqual(loaded_base.molecule, sidecar_base.molecule)
            self.assertEqual(loaded_base.nucleotide, sidecar_base.nucleotide)
            self.assertEqual(
                coordinate_map.get(loaded_base.up),
                sidecar_base.up,
            )
            self.assertEqual(
                coordinate_map.get(loaded_base.down),
                sidecar_base.down,
            )
            self.assertEqual(
                coordinate_map.get(loaded_base.across),
                sidecar_base.across,
            )

    def test_rna_dna_binary_loads_with_rna_record_layout(self) -> None:
        if not SAMPLE_RNA_DNA_PATH.exists() or not SAMPLE_RNA_DNAJSON_PATH.exists():
            self.skipTest("RNA .dna/.dnajson fixtures are not available")
        records = parse_dna_records(SAMPLE_RNA_DNA_PATH)
        project = load_dna_project(SAMPLE_RNA_DNA_PATH)
        sidecar = load_json_project(SAMPLE_RNA_DNAJSON_PATH)
        self.assertEqual(len(records), 1680)
        self.assertEqual(records[0].nucleotide, "U")
        self.assertEqual(records[0].molecule, "RNA")
        self.assertEqual(project.summary(), sidecar.summary())
        self.assertEqual(project.metadata.get("dna_sidecar_path"), str(SAMPLE_RNA_DNAJSON_PATH))


if __name__ == "__main__":
    unittest.main()
