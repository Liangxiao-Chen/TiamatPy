from __future__ import annotations

import argparse
from pathlib import Path
import sys

try:
    from .io_ascii import load_ascii_dat
    from .io_dna import load_dna_project
    from .io_json import load_json_project, save_json_project
    from .io_pdb import load_pdb_sugar_center_project
    from .ncbi import fetch_fasta, search_genomes
    from .viewer import Camera, launch, load_project, render_svg
except ImportError:
    from tiamat_py.io_ascii import load_ascii_dat
    from tiamat_py.io_dna import load_dna_project
    from tiamat_py.io_json import load_json_project, save_json_project
    from tiamat_py.io_pdb import load_pdb_sugar_center_project
    from tiamat_py.ncbi import fetch_fasta, search_genomes
    from tiamat_py.viewer import Camera, launch, load_project, render_svg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pragmatic Python port of core Tiamat workflows.")
    subparsers = parser.add_subparsers(dest="command")

    gui = subparsers.add_parser("gui", help="open the Tk viewer")
    gui.add_argument("input", nargs="?", help="optional .dna, .dnajson, or .pdb file")

    inspect = subparsers.add_parser("inspect", help="print a project summary")
    inspect.add_argument("input", help=".dna, .dnajson, or .pdb file")

    convert = subparsers.add_parser("convert", help="convert .dat or .pdb to .dnajson")
    convert.add_argument("input", help="input .dat or .pdb file")
    convert.add_argument("output", help="output .dnajson file")

    export = subparsers.add_parser("export-sequences", help="write strand sequences to text")
    export.add_argument("input", help=".dna, .dnajson, or .pdb file")
    export.add_argument("output", help="output .txt file")
    export.add_argument("--generate", action="store_true", help="generate sequences before export")

    svg = subparsers.add_parser("export-svg", help="render a static SVG snapshot")
    svg.add_argument("input", help=".dna, .dnajson, or .pdb file")
    svg.add_argument("output", help="output .svg file")

    search = subparsers.add_parser("search-genome", help="search NCBI genome records")
    search.add_argument("term", help="NCBI genome query")
    search.add_argument("--limit", type=int, default=5, help="maximum matches to return")

    fasta = subparsers.add_parser("fetch-fasta", help="download FASTA for a genome id")
    fasta.add_argument("genome_id", help="NCBI genome id")
    fasta.add_argument("--output", help="optional output .fasta file")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command in (None, "gui"):
        project = load_project(args.input) if getattr(args, "input", None) else None
        launch(project)
        return 0

    if args.command == "inspect":
        project = load_project(args.input)
        summary = project.summary()
        print(f"bases={summary['bases']}")
        print(f"strands={summary['strands']}")
        print(f"paired_bases={summary['paired_bases']}")
        print(f"assigned_bases={summary['assigned_bases']}")
        return 0

    if args.command == "convert":
        input_path = Path(args.input)
        if input_path.suffix.lower() == ".dna":
            project = load_dna_project(input_path)
        elif input_path.suffix.lower() == ".pdb":
            project = load_pdb_sugar_center_project(input_path)
        else:
            project = load_ascii_dat(input_path)
        save_json_project(project, args.output)
        print(args.output)
        return 0

    if args.command == "export-sequences":
        project = load_project(args.input)
        if args.generate:
            project.generate_sequences()
        Path(args.output).write_text(project.export_sequences_text(), encoding="utf-8")
        print(args.output)
        return 0

    if args.command == "export-svg":
        project = load_project(args.input)
        svg = render_svg(project, Camera())
        Path(args.output).write_text(svg, encoding="utf-8")
        print(args.output)
        return 0

    if args.command == "search-genome":
        for match in search_genomes(args.term, limit=args.limit):
            print(f"{match.genome_id}\t{match.title}")
        return 0

    if args.command == "fetch-fasta":
        fasta_text = fetch_fasta(args.genome_id)
        if args.output:
            Path(args.output).write_text(fasta_text, encoding="utf-8")
            print(args.output)
        else:
            sys.stdout.write(fasta_text)
        return 0

    parser.error(f"unsupported command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
