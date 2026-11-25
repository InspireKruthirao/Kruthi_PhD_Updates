"""Helpers shared by the Supacow gene catalogue scripts."""

import gzip
import lzma
import os
import subprocess


def grouped_id(prefix, level, ix):
    """SC.ORF.000_000_001 style ID (Shanghai Dogs / Urban Soil convention)"""
    n = f"{ix:09d}"
    return f"{prefix}.{level}.{n[:3]}_{n[3:6]}_{n[6:]}"


def id_number(identifier):
    """SC.100NT.000_012_345 -> 12345"""
    return int(identifier.rsplit(".", 1)[1].replace("_", ""))


def open_text(path):
    if str(path).endswith(".xz"):
        return lzma.open(path, "rt")
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fasta_iter(path):
    """(full header, sequence) for each record"""
    header, seq = None, []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


class AtomicWriter:
    """Write to <path>.part (xz-compressed with threads if path ends in .xz,
    gzip if .gz) and rename to <path> only when closed without error."""

    def __init__(self, path, threads=8):
        self.path, self.part = str(path), str(path) + ".part"
        if os.path.exists(self.part):
            os.remove(self.part)
        self.proc = None
        if self.path.endswith(".xz"):
            self.raw = open(self.part, "wb")
            self.proc = subprocess.Popen(["xz", f"-T{threads}", "-c"], stdin=subprocess.PIPE, stdout=self.raw)
            self.f = self.proc.stdin
            self.binary = True
        elif self.path.endswith(".gz"):
            self.f, self.binary = gzip.open(self.part, "wb"), True
        else:
            self.f, self.binary = open(self.part, "w"), False

    def write(self, text):
        self.f.write(text.encode() if self.binary else text)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.f.close()
        if self.proc is not None:
            status = self.proc.wait()
            self.raw.close()
            if status != 0 and exc_type is None:
                raise RuntimeError(f"xz failed writing {self.path}")
        if exc_type is None:
            os.replace(self.part, self.path)
        return False


def check(condition, message):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + message)
    print("  ok: " + message, flush=True)
