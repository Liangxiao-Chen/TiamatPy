from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class GenomeMatch:
    genome_id: str
    title: str


def search_genomes(term: str, limit: int = 5) -> list[GenomeMatch]:
    params = urlencode(
        {
            "tool": "TiamatPy",
            "db": "genome",
            "term": term,
            "retmax": limit,
        }
    )
    root = ET.fromstring(_fetch(f"{EUTILS_BASE}/esearch.fcgi?{params}"))
    genome_ids = [element.text for element in root.findall("./IdList/Id") if element.text]
    return [GenomeMatch(genome_id=genome_id, title=_fetch_fasta_header(genome_id)) for genome_id in genome_ids]


def fetch_fasta(genome_id: str) -> str:
    params = urlencode(
        {
            "tool": "TiamatPy",
            "db": "genome",
            "retmode": "text",
            "rettype": "fasta",
            "id": genome_id,
        }
    )
    return _fetch(f"{EUTILS_BASE}/efetch.fcgi?{params}")


def _fetch_fasta_header(genome_id: str) -> str:
    params = urlencode(
        {
            "tool": "TiamatPy",
            "db": "genome",
            "retmode": "text",
            "rettype": "fasta",
            "seq_stop": 1,
            "id": genome_id,
        }
    )
    text = _fetch(f"{EUTILS_BASE}/efetch.fcgi?{params}")
    first_line = text.splitlines()[0] if text else genome_id
    return first_line.lstrip(">")


def _fetch(url: str) -> str:
    request = Request(url, headers={"User-Agent": "TiamatPy/0.1"})
    with urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8", errors="ignore")

