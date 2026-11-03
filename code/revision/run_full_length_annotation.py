"""Run the archived EBI client on corrected full-length protein sequences.

Usage: python run_full_length_annotation.py LEGACY_CLIENT FULL_FASTA OUTPUT_DIR
The client is loaded without invoking its original main function. No historical
annotation output is overwritten. The EBI endpoint and returned job IDs are saved.
"""
import importlib.util
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    legacy, fasta, output = map(Path, sys.argv[1:])
    output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("archived_ebi_client", legacy)
    client = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(client)
    client.FASTA = str(fasta)
    client.OUTDIR = str(output)
    client.STATE = str(output / "iprscan_jobs_state.json")
    original_http = client.http
    def guarded_http(*args, **kwargs):
        try:
            return original_http(*args, **kwargs)
        except RuntimeError as error:
            if "blocked" in str(error).lower():
                raise SystemExit("EBI submission blocked by provider fair-use policy; no further submissions.") from error
            raise
    client.http = guarded_http
    client.main()


if __name__ == "__main__":
    main()
