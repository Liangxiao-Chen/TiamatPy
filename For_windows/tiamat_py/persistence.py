"""Failure-safe replacement of an individual document on the same filesystem."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile
import shutil


def atomic_write_text(path: str | Path, text: str) -> None:
    destination = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         dir=destination.parent, prefix=f".{destination.name}.",
                                         suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        if destination.exists():
            temporary.chmod(destination.stat().st_mode & 0o777)
        os.replace(temporary, destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def atomic_write_files(contents: dict[Path, str]) -> None:
    """Stage all outputs, then replace; roll back a partial ordinary I/O failure.

    Filesystems cannot atomically replace two separate paths across power loss.
    A failed rollback deliberately leaves the recovery directory on disk.
    """
    paths = list(contents)
    staging = Path(tempfile.mkdtemp(prefix=".tiamat-export-", dir=paths[0].parent))
    replaced = []
    backups = {}
    try:
        for index, path in enumerate(paths):
            if path.exists():
                backup = staging / f"{index}.backup"
                shutil.copy2(path, backup)
                backups[path] = backup
            atomic_write_text(staging / str(index), contents[path])
        for index, path in enumerate(paths):
            os.replace(staging / str(index), path)
            replaced.append(path)
    except BaseException:
        try:
            for path in reversed(replaced):
                if path in backups:
                    os.replace(backups[path], path)
                else:
                    path.unlink(missing_ok=True)
        except OSError as recovery_error:
            raise OSError(f"Export rollback failed. Original files are preserved in {staging}") from recovery_error
        shutil.rmtree(staging)
        raise
    shutil.rmtree(staging)
