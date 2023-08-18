"""Mock test URL download."""

import io
import os

from gencdna.file_api.parse_exons_introns import write_exons_introns


def test_app_exons_and_introns_fasta_from_url(
    app_gene_html,
    ref_app_exons_fasta,
    snapshot,
):
    with io.StringIO() as app_gene_exons_and_introns_fasta:
        snapshot.snapshot_dir = os.path.dirname(ref_app_exons_fasta)
        write_exons_introns(
            ucsc_html=app_gene_html,
            exons_fasta=app_gene_exons_and_introns_fasta,
            introns_fasta='',
        )
        snapshot.assert_match(
            app_gene_exons_and_introns_fasta.getvalue(),
            ref_app_exons_fasta,
        )
