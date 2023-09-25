Next steps
==========

Finish finding exon-exon joins for all samples - last step in
:ref:`wgs-bwa-pipeline`.

Summarize results
-----------------

Starting from ``<mapped_exons_joins_csv>`` files and create a summary of all
exon-exon joins found across the sample.


Human genome insertion sites
----------------------------

For each read where we found exon-exon joins, which is two exons on a read
without a gap, check if this is already present in a reference human genome. 
Generally, we are interested in somatic and not germline exon-exon joins so 
these need to be filtered out.

Besides creating a custom code for this, it's worth checking other tools such
as `Sniffles <https://github.com/fritzsedlazeck/Sniffles>`_. 
