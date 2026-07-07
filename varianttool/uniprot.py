from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import quote
import re

import requests

FASTA_BASE = "https://rest.uniprot.org/uniprotkb/{source}.fasta"
STREAM_BASE = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query={query}"


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str
    url: str

    @property
    def gene_symbol(self) -> str:
        m = re.search(r"\bGN=([^\s]+)", self.header)
        return m.group(1) if m else ""


def _is_accession(source: str) -> bool:
    return (" " not in source) and ("AND" not in source) and (":" not in source)


def _fetch_text(source: str, timeout: int = 45) -> tuple[str, str]:
    if _is_accession(source):
        url = FASTA_BASE.format(source=source)
    else:
        url = STREAM_BASE.format(query=quote(source, safe=""))
    r = requests.get(url, timeout=timeout, headers={"User-Agent": "varianttool/0.1.0"})
    r.raise_for_status()
    return r.text, url


def fetch_fasta_records(source: str, timeout: int = 45) -> list[FastaRecord]:
    text, url = _fetch_text(source, timeout=timeout)
    records: list[FastaRecord] = []
    header = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(FastaRecord(header=header, sequence="".join(seq_parts).upper(), url=url))
            header = line[1:].strip()
            seq_parts = []
        else:
            if header is not None:
                seq_parts.append(line.strip())
    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(seq_parts).upper(), url=url))
    if not records:
        raise ValueError(f"Could not parse FASTA from {source!r}")
    return records


def fetch_fasta(source: str, timeout: int = 45) -> tuple[str, str, str]:
    rec = fetch_fasta_records(source, timeout=timeout)[0]
    return rec.header, rec.sequence, rec.url
