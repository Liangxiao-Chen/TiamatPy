"""Portable regression cases for the eight 2026-10 audit findings."""
from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import Mock, patch

from tiamat_py.dna_archive import decode_archive, encode_archive
from tiamat_py.geometry import add, cross, dot, norm, scale
from tiamat_py.io_dna import load_dna_project, parse_dna_records
from tiamat_py.io_json import load_json_project, save_json_project
from tiamat_py.io_oxdna import BACKBONE_OFFSETS, NM_PER_OXDNA_UNIT, build_oxdna_configuration, export_oxdna_files
from tiamat_py.io_pdb import load_pdb_sugar_center_project, save_minimal_dna_project
from tiamat_py.model import Base, TiamatProject, basepair_state
from tiamat_py.persistence import atomic_write_files, atomic_write_text
from tiamat_py.viewer import (TiamatViewer, _apply_rotation_to_project, _blank_custom_constraints_values,
    _build_created_structure, _default_constraints_values, _snapshot_history_state, _restore_history_state)


def chain(count=12, molecule="DNA", circular=False):
    return TiamatProject({i: Base(i, str(i), i * .6, 0, 0, molecule=molecule,
        up=(i-1) % count if i > 0 or circular else None,
        down=(i+1) % count if i < count-1 or circular else None) for i in range(count)})


def helix(count=20, molecule="DNA-DNA B", kind="Helix"):
    return TiamatProject(_build_created_structure((0, 0, 0), (1, 0, 0), count, 0, molecule,
                        kind, _default_constraints_values(), 0))


def headless_viewer(project):
    viewer = TiamatViewer.__new__(TiamatViewer)
    viewer.project = project
    viewer.root = Mock()
    viewer.document_path = None
    viewer.constraints_values = _default_constraints_values()
    viewer.custom_constraints_values = _blank_custom_constraints_values()
    viewer._saved_document = viewer._document_content()
    viewer.selected_indices = set()
    viewer.primary_selected = None
    viewer.undo_stack = []
    viewer.redo_stack = []
    for method in ("_refresh_all", "_remember_recent_file", "_clear_pending_free_strand", "_reset_viewports"):
        setattr(viewer, method, Mock())
    return viewer


def pdb_residue(number=1, name="DA", chain_id="A", center=(-.6, 0, 0), right=False,
                shift=(0, 0, 0), alt=" ", occupancy=1, sugar_only=False, wobble=False):
    # Synthetic planar bases with two 0.29-nm Watson-Crick contacts.
    atoms = {"C1'": center, "C2'": add(center, (.05, 0, 0)), "C3'": add(center, (0, .05, 0)),
             "C4'": add(center, (-.05, 0, 0)), "O4'": add(center, (0, -.05, 0))}
    if not sugar_only:
        if right:
            atoms.update({"N1": (.6, -.14, 0), "C2": (.4, -.1, 0), "N3": (.29, 0, 0),
                          "C4": (.45, .14, 0), "C5": (.55, .3, 0), "C6": (.6, 0, 0), "O4": (.29, .3, 0)})
        else:
            atoms.update({"N1": (0, 0, 0), "C2": (0, .14, 0), "N3": (-.1, .2, 0),
                          "C4": (-.24, .14, 0), "C5": (-.24, 0, 0), "C6": (-.12, -.07, 0), "N6": (0, .3, 0)})
        if wobble:
            if right:
                atoms["O2"] = (.29, 0, 0)
                atoms["N3"] = (.29, .3, 0)
            else:
                atoms["O6"] = atoms.pop("N6")
    lines = []
    for serial, (atom, point) in enumerate(atoms.items(), 1):
        x, y, z = scale(add(point, shift), 10)
        lines.append(f"ATOM  {serial:5d} {atom:>4s}{alt}{name:>3s} {chain_id}{number:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{0:6.2f}\n")
    return "".join(lines)


class PersistenceTests(unittest.TestCase):
    def test_zero_id_survives_native_json_and_history(self):
        project = helix(4)
        state = _snapshot_history_state(project, {0}, 0)
        restored, selected, primary = _restore_history_state(state)
        self.assertEqual(restored.to_dict(), project.to_dict())
        self.assertEqual((selected, primary), ({0}, 0))
        self.assertEqual(restored.bases[1].across, 0)
        self.assertEqual(restored.bases[2].up, 0)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "zero.dnajson"
            save_json_project(project, path)
            loaded = load_json_project(path)
            self.assertEqual(loaded.bases[1].across, 0)
            self.assertEqual(loaded.bases[2].up, 0)

    def test_circular_zero_id_round_trip(self):
        project = chain(6, circular=True)
        restored = TiamatProject.from_dict(project.to_dict())
        self.assertEqual(restored.bases[5].down, 0)
        self.assertEqual(restored.strands[0].base_indices, list(range(6)))

    def test_save_as_reopen_rename_uses_current_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            original, copied, renamed = [Path(tmp) / name for name in ("original.dnajson", "copy.dnajson", "renamed.dnajson")]
            viewer = headless_viewer(chain())
            self.assertTrue(viewer._write_project_json(original))
            before = original.read_bytes()
            viewer.project.bases[0].x = 50
            self.assertTrue(viewer._write_project_json(copied))
            self.assertEqual(json.loads(copied.read_text())["metadata"]["source_path"], str(copied.resolve()))
            copied.rename(renamed)
            loaded = load_json_project(renamed)
            self.assertEqual(loaded.metadata["source_path"], str(renamed.resolve()))
            reopened = headless_viewer(loaded)
            reopened.document_path = renamed
            reopened.project.bases[0].x = 60
            self.assertTrue(reopened.save_project())
            self.assertEqual(original.read_bytes(), before)
            self.assertEqual(load_json_project(renamed).bases[0].x, 60)

    def test_undo_cannot_redirect_save_path(self):
        viewer = headless_viewer(chain())
        snapshot = _snapshot_history_state(viewer.project, set(), None)
        viewer.document_path = Path("/tmp/new-name.dnajson")
        viewer._restore_history_snapshot(snapshot)
        self.assertEqual(viewer._project_save_path(), Path("/tmp/new-name.dnajson"))

    def test_selection_is_not_unsaved_edit(self):
        viewer = headless_viewer(chain())
        viewer.selected_indices = {0, 1}
        with patch("tiamat_py.viewer.messagebox.askyesnocancel") as prompt:
            self.assertTrue(viewer._confirm_discard_changes())
            prompt.assert_not_called()

    def test_cancel_prevents_new_and_close(self):
        viewer = headless_viewer(chain())
        project = viewer.project
        project.bases[0].x = 12
        with patch("tiamat_py.viewer.messagebox.askyesnocancel", return_value=None):
            viewer._new_project()
            viewer._close_window()
        self.assertIs(viewer.project, project)
        viewer.root.destroy.assert_not_called()

    def test_save_dialog_cancel_prevents_close(self):
        viewer = headless_viewer(chain())
        viewer.project.bases[0].x = 12
        with patch("tiamat_py.viewer.messagebox.askyesnocancel", return_value=True), \
             patch("tiamat_py.viewer.filedialog.asksaveasfilename", return_value=""):
            viewer._close_window()
        viewer.root.destroy.assert_not_called()

    def test_failed_save_keeps_dirty_path_and_document(self):
        viewer = headless_viewer(chain())
        viewer.project.bases[0].x = 12
        with patch("tiamat_py.viewer.save_json_project", side_effect=OSError("disk full")), \
             patch("tiamat_py.viewer.messagebox.showerror"):
            self.assertFalse(viewer._write_project_json("new.dnajson"))
        self.assertIsNone(viewer.document_path)
        self.assertNotEqual(viewer._saved_document, viewer._document_content())

    def test_discard_allows_new(self):
        viewer = headless_viewer(chain())
        viewer.project.bases[0].x = 12
        with patch("tiamat_py.viewer.messagebox.askyesnocancel", return_value=False):
            viewer._new_project()
        self.assertIsNone(viewer.project)

    def test_open_failure_keeps_current_document(self):
        viewer = headless_viewer(chain())
        before = viewer.project
        with self.assertRaises(FileNotFoundError):
            viewer._open_project_path("/no-such-tiamat-file.dnajson")
        self.assertIs(viewer.project, before)

    def test_open_same_file_after_save_uses_new_contents(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"current.dnajson"
            viewer = headless_viewer(chain())
            viewer._write_project_json(path)
            viewer.project.bases[0].x = 42
            with patch("tiamat_py.viewer.messagebox.askyesnocancel", return_value=True):
                viewer._open_project_path(path)
            self.assertEqual(viewer.project.bases[0].x, 42)

    def test_open_cancel_does_not_read_or_replace_document(self):
        viewer = headless_viewer(chain())
        viewer.project.bases[0].x = 42
        with patch("tiamat_py.viewer.messagebox.askyesnocancel", return_value=None), \
             patch("tiamat_py.viewer.load_project") as loader:
            viewer._open_project_path("/not-read.dnajson")
        loader.assert_not_called()
        self.assertEqual(viewer.project.bases[0].x, 42)

    def test_atomic_save_keeps_original_and_cleans_temp_on_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "document.dnajson"
            path.write_text("original")
            with patch("tiamat_py.persistence.os.replace", side_effect=OSError("failure")):
                with self.assertRaises(OSError):
                    atomic_write_text(path, "replacement")
            self.assertEqual(path.read_text(), "original")
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_two_file_export_rolls_back_second_replace_failure(self):
        import os
        original_replace = os.replace
        with tempfile.TemporaryDirectory() as tmp:
            top, conf = Path(tmp)/"data.top", Path(tmp)/"data.oxdna"
            top.write_text("old top")
            conf.write_text("old conf")
            def replace(source, destination):
                if Path(destination) == conf:
                    raise OSError("second file failed")
                return original_replace(source, destination)
            with patch("tiamat_py.persistence.os.replace", side_effect=replace):
                with self.assertRaises(OSError):
                    atomic_write_files({top: "new top", conf: "new conf"})
            self.assertEqual(top.read_text(), "old top")
            self.assertEqual(conf.read_text(), "old conf")


class GraphTests(unittest.TestCase):
    def test_missing_neighbor_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "missing base 9"):
            TiamatProject({0: Base(0, "0", 0, 0, 0, up=9)})

    def test_nonreciprocal_links_are_rejected(self):
        for link in ("up", "down", "across"):
            with self.subTest(link=link), self.assertRaisesRegex(ValueError, "not reciprocal"):
                TiamatProject({0: Base(0, "0", 0, 0, 0, **{link: 1}), 1: Base(1, "1", 1, 0, 0)})

    def test_self_links_and_nonfinite_coordinates_rejected(self):
        for kwargs in ({"up": 0}, {"across": 0}, {"x": float("nan")}, {"z": float("inf")}):
            values = dict(index=0, object_id="0", x=0, y=0, z=0)
            values.update(kwargs)
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                TiamatProject({0: Base(**values)})

    def test_duplicate_ids_rejected_in_both_json_formats(self):
        native = chain(2).to_dict()
        native["bases"].append(dict(native["bases"][0]))
        external = {"bases": [{"id": 0, "position": [0, 0, 0]}, {"id": "0", "position": [1, 0, 0]}]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"bad.dnajson"
            for data in (native, external):
                path.write_text(json.dumps(data))
                with self.assertRaisesRegex(ValueError, "[Dd]uplicate"):
                    load_json_project(path)

    def test_mutated_graph_export_fails_instead_of_looping(self):
        project = helix(4)
        project.bases[0].up = 0
        start = time.monotonic()
        with self.assertRaises(ValueError):
            export_oxdna_files(project, "should-not-exist.top", default_base="R")
        self.assertLess(time.monotonic()-start, 1)

    def test_invalid_molecule_letter_edit_is_atomic(self):
        project = helix(2)
        before = project.to_dict()
        with self.assertRaises(ValueError):
            project.set_nucleotide(0, "U")
        self.assertEqual(project.to_dict(), before)


class ArchiveTests(unittest.TestCase):
    def test_repository_demo_matches_json_without_neighbors(self):
        root = Path(__file__).resolve().parents[1]
        reference = load_json_project(root/"demo/test.dnajson")
        with tempfile.TemporaryDirectory() as tmp:
            dna = Path(tmp)/"moved.dna"
            dna.write_bytes((root/"demo/test.dna").read_bytes())
            parsed = load_dna_project(dna)
            save_json_project(chain(3), dna.with_suffix(".dnajson"))
            again = load_dna_project(dna)
            self.assertEqual(parsed.to_dict(), again.to_dict())
        self.assertEqual(parsed.summary(), reference.summary())
        for index, base in parsed.bases.items():
            target = reference.bases[index]
            self.assertLess(math.dist(base.position, target.position), .01)
            self.assertEqual((base.molecule, base.nucleotide, base.up, base.down, base.across),
                             (target.molecule, target.nucleotide, target.up, target.down, target.across))

    def test_generic_rna_at_origin_and_circular_topology(self):
        project = chain(8, molecule="RNA", circular=True)
        decoded, records = decode_archive(encode_archive(project))
        self.assertEqual(len(records), 8)
        for i in project.bases:
            self.assertEqual(decoded.bases[i].to_dict(), project.bases[i].to_dict())

    def test_schema5_roundtrip_dna_rna_and_hybrid(self):
        for molecule in ("DNA-DNA B", "RNA-RNA A", "DNA-RNA A"):
            with self.subTest(molecule=molecule):
                project = helix(10, molecule)
                project.generate_sequences(seed=3)
                decoded, _ = decode_archive(encode_archive(project))
                for i, base in project.bases.items():
                    target = decoded.bases[i]
                    self.assertEqual((base.position, base.molecule, base.nucleotide, base.up, base.down, base.across),
                                     (target.position, target.molecule, target.nucleotide, target.up, target.down, target.across))

    def test_deep_archive_uses_stack_not_python_recursion(self):
        project = chain(4000)
        decoded, _ = decode_archive(encode_archive(project))
        self.assertEqual(decoded.bases[3999].up, 3998)

    def test_extended_object_references_above_32767(self):
        project = chain(33000)
        decoded, _ = decode_archive(encode_archive(project))
        self.assertEqual(decoded.bases[32999].up, 32998)

    def test_truncation_and_wrong_schema_rejected(self):
        data = encode_archive(helix(3))
        for broken in (data[:-1], data[:250], data[:2]+b"\x06\x00"+data[4:]):
            with self.subTest(length=len(broken)), self.assertRaises(ValueError):
                decode_archive(broken)

    def test_pdb_archive_writer_needs_no_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = helix(4)
            path, sidecar = save_minimal_dna_project(project, Path(tmp)/"made.dna", include_sidecar=False)
            self.assertIsNone(sidecar)
            decoded = load_dna_project(path)
            self.assertEqual(decoded.summary(), project.summary())
            self.assertEqual(len(parse_dna_records(path)), 8)


class PdbTests(unittest.TestCase):
    def load_text(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"synthetic.pdb"
            path.write_text(text)
            return load_pdb_sugar_center_project(path)

    def test_atom_geometry_pairs_close_complementary_bases(self):
        project = self.load_text(pdb_residue() + pdb_residue(name="DT", chain_id="B", center=(.9, 0, 0), right=True))
        self.assertEqual(project.bases[1].across, 2)
        self.assertEqual(project.bases[2].across, 1)
        self.assertAlmostEqual(project.bases[1].x, -.6)

    def test_far_or_nonplanar_complements_are_not_paired(self):
        for shift in ((100, 0, 0), (0, 0, .4)):
            project = self.load_text(pdb_residue() + pdb_residue(name="DT", chain_id="B", right=True, shift=shift))
            self.assertIsNone(project.bases[1].across)

    def test_gu_wobble_uses_rna_contacts(self):
        project = self.load_text(pdb_residue(name="G", wobble=True)
            + pdb_residue(name="U", chain_id="B", center=(.9, 0, 0), right=True, wobble=True))
        self.assertEqual(project.bases[1].across, 2)
        self.assertEqual(basepair_state(project.bases[1], project.bases[2]), "wobble")

    def test_ambiguous_geometric_partners_are_left_unpaired(self):
        project = self.load_text(pdb_residue()
            + pdb_residue(name="DT", chain_id="B", center=(.9, 0, 0), right=True)
            + pdb_residue(name="DT", chain_id="C", center=(.9, 0, 0), right=True))
        self.assertEqual(project.summary()["paired_bases"], 0)

    def test_more_than_two_chains_and_hairpin_pairs(self):
        text = pdb_residue() + pdb_residue(number=8, name="DT", right=True, center=(.9, 0, 0))
        text += pdb_residue(chain_id="C", shift=(0, 3, 0))
        text += pdb_residue(name="DT", chain_id="D", right=True, center=(.9, 0, 0), shift=(0, 3, 0))
        project = self.load_text(text)
        self.assertEqual(project.summary()["paired_bases"], 4)

    def test_sugar_only_file_does_not_invent_pairs(self):
        project = self.load_text(pdb_residue(sugar_only=True) + pdb_residue(name="DT", chain_id="B", sugar_only=True))
        self.assertEqual(project.summary()["paired_bases"], 0)

    def test_ter_prevents_backbone_connection(self):
        left = pdb_residue(center=(0, 0, 0), sugar_only=True)
        right = pdb_residue(number=2, center=(.6, 0, 0), sugar_only=True)
        self.assertEqual(len(self.load_text(left + right).strands), 1)
        self.assertEqual(len(self.load_text(left + "TER\n" + right).strands), 2)

    def test_first_model_only(self):
        text = "MODEL        1\n" + pdb_residue() + "ENDMDL\nMODEL        2\n" + pdb_residue(shift=(100, 0, 0)) + "ENDMDL\n"
        self.assertAlmostEqual(self.load_text(text).bases[1].x, -.6)

    def test_altloc_uses_one_highest_occupancy_conformer(self):
        text = pdb_residue(alt="A", occupancy=.8) + pdb_residue(alt="B", occupancy=.2, shift=(10, 0, 0))
        self.assertAlmostEqual(self.load_text(text).bases[1].x, -.6)


class SequenceTests(unittest.TestCase):
    def assert_constraints(self, project, unique=10, repeat=6):
        for strand in project.strands:
            seq = strand.sequence(project).replace("U", "T")
            circular = project.bases[strand.base_indices[-1]].down == strand.base_indices[0]
            if circular:
                seq += (seq * (max(unique, repeat)//len(seq)+1))[:max(unique, repeat)-1]
            for start in range(len(seq)-unique+1):
                self.assertEqual(set(seq[start:start+unique]), set("ATCG"))
            for start in range(len(seq)-repeat+1):
                self.assertGreater(len(set(seq[start:start+repeat])), 1)
        for base in project.bases.values():
            if base.across is not None:
                self.assertEqual(basepair_state(base, project.bases[base.across]), "canonical")
        fraction = sum(b.nucleotide in "GC" for b in project.bases.values())/len(project.bases)
        self.assertTrue(.45 <= fraction <= .55)

    def test_default_is_not_period_four(self):
        sequences = []
        for seed in range(5):
            project = chain(64)
            project.generate_sequences(seed=seed)
            self.assert_constraints(project)
            seq = project.strands[0].sequence(project)
            self.assertTrue(any(seq[i] != seq[i+4] for i in range(len(seq)-4)))
            sequences.append(seq)
        self.assertEqual(len(set(sequences)), 5)

    def test_seed_is_reproducible(self):
        first, second = helix(40), helix(40)
        first.generate_sequences(seed=21)
        second.generate_sequences(seed=21)
        self.assertEqual(first.export_sequences_text(), second.export_sequences_text())

    def test_dna_rna_hybrid_and_circular_constraints(self):
        for project in (helix(50), helix(50, "RNA-RNA A"), helix(50, "DNA-RNA A"), chain(32, circular=True)):
            project.generate_sequences(seed=18)
            self.assert_constraints(project)

    def test_timeout_keeps_original_sequence(self):
        project = chain(5, circular=True)
        for base in project.bases.values():
            base.nucleotide = "A"
        before = project.to_json()
        with self.assertRaises(TimeoutError):
            project.generate_sequences(unique_sequence_limit=4, gc_tolerance=.5, timeout=.02)
        self.assertEqual(project.to_json(), before)

    def test_large_input_has_bounded_runtime(self):
        project = helix(1000)
        project.generate_sequences(seed=4, timeout=5)
        self.assert_constraints(project)


class OxDnaTests(unittest.TestCase):
    def assert_vector_close(self, left, right):
        self.assertLess(math.dist(left, right), 1e-9)

    def test_frames_orthonormal_and_recover_original_backbones(self):
        for molecule, model in (("DNA-DNA B", "DNA2"), ("RNA-RNA A", "RNA")):
            project = helix(20, molecule)
            configs = build_oxdna_configuration(project)
            for index, (center, a1, a3) in configs.items():
                self.assertAlmostEqual(norm(a1), 1)
                self.assertAlmostEqual(norm(a3), 1)
                self.assertAlmostEqual(dot(a1, a3), 0)
                a2 = cross(a3, a1)
                offsets = BACKBONE_OFFSETS[model]
                backbone = add(center, add(add(scale(a1, offsets[0]), scale(a2, offsets[1])), scale(a3, offsets[2])))
                self.assert_vector_close(scale(backbone, NM_PER_OXDNA_UNIT), project.bases[index].position)

    def test_rigid_transform_equivariance_with_unpaired_tails(self):
        for molecule, kind in (("DNA-DNA B", "Helix"), ("RNA-RNA A", "Helix"), ("DNA-DNA B", "Strand")):
            project = helix(12, molecule, kind)
            if kind == "Helix":
                project.bases[0].across = project.bases[1].across = None
            before = build_oxdna_configuration(project)
            def rotate(v): return (-v[1], v[2], -v[0])
            translation = (20, 30, 40)
            for base in project.bases.values():
                base.x, base.y, base.z = add(rotate(base.position), translation)
                base.plane_u, base.plane_v = rotate(base.plane_u), rotate(base.plane_v)
            after = build_oxdna_configuration(project)
            for index, (center, a1, a3) in before.items():
                self.assert_vector_close(after[index][0], add(rotate(center), scale(translation, 1/NM_PER_OXDNA_UNIT)))
                self.assert_vector_close(after[index][1], rotate(a1))
                self.assert_vector_close(after[index][2], rotate(a3))

    def test_unpaired_tail_without_stored_frames_transports_from_duplex(self):
        project = helix(12)
        project.bases[0].across = project.bases[1].across = None
        for base in project.bases.values():
            base.plane_u = base.plane_v = None
        before = build_oxdna_configuration(project)
        for base in project.bases.values():
            base.x += 500
        after = build_oxdna_configuration(project)
        for index in before:
            self.assert_vector_close(before[index][1], after[index][1])

    def test_underdetermined_straight_strand_reports_error(self):
        with self.assertRaisesRegex(ValueError, "no defined base plane"):
            build_oxdna_configuration(chain())

    def test_drag_rotation_uses_original_frames_not_cumulative_angle(self):
        project = helix(4)
        positions = {i: b.position for i, b in project.bases.items()}
        frames = {i: (b.plane_u, b.plane_v) for i, b in project.bases.items()}
        rotation = ((0., -1., 0.), (1., 0., 0.), (0., 0., 1.))
        _apply_rotation_to_project(project, positions, (0, 0, 0), rotation, original_frames=frames)
        once = project.to_dict()
        _apply_rotation_to_project(project, positions, (0, 0, 0), rotation, original_frames=frames)
        self.assertEqual(once, project.to_dict())

    def test_random_default_keeps_pairs_and_classic_direction(self):
        project = helix(20)
        with tempfile.TemporaryDirectory() as tmp:
            top, conf = export_oxdna_files(project, Path(tmp)/"helix", default_base="R")
            rows = [line.split() for line in top.read_text().splitlines()[1:]]
            order = [i for strand in project.strands for i in reversed(strand.base_indices)]
            lookup = {index: row for index, row in zip(order, rows)}
            for index, base in project.bases.items():
                val = lookup[index][1]
                self.assertEqual(lookup[base.across][1], {"A":"T", "T":"A", "G":"C", "C":"G"}[val])
                row = lookup[index]
                self.assertEqual(order[int(row[2])] if row[2] != "-1" else None, base.down)
                self.assertEqual(order[int(row[3])] if row[3] != "-1" else None, base.up)
            self.assertTrue(all(base.nucleotide is None for base in project.bases.values()))
            self.assertEqual(len(conf.read_text().splitlines()), 43)

    def test_large_structure_box_contains_all_centers(self):
        project = helix(4)
        for base in project.bases.values():
            base.x += 1000
        with tempfile.TemporaryDirectory() as tmp:
            _, conf = export_oxdna_files(project, Path(tmp)/"large", default_base="R")
            lines = conf.read_text().splitlines()
            box = float(lines[1].split()[2])
            for line in lines[3:]:
                self.assertLess(max(abs(float(v)) for v in line.split()[:3]), box/2)


if __name__ == "__main__":
    unittest.main()
