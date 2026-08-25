# Bugs fixed relative to upstream RIsearch 1.2

Defects in upstream [RTH-tools/risearch](https://github.com/RTH-tools/risearch) that
this fork fixes. Each was reproduced against the original C and has a regression test.

* Fixed infinite loop with FASTA input files where the last line is a header
* A record with no sequence, or one of only gap characters, will overflow
* Alignment with no complementary pair would access invalid memory
* Corrupt bytes in a sequence were dropped without a word, so a damaged record was
  read as a shorter one that still looked valid, and scored as if it were whole
* Fixed isalpha (undefined behavior for bytes above 127, since a plain char is signed)
* Switched isalpha with a table lookup, improved correctness and FASTA reading speeds
* Leaked memory when sequence was rejected `-Q "ACGU@ACGU"`

Behavioral changes:

* Sequence with a control byte, or a byte of value 127 or above, is rejected instead
  of omitted
* A record that encodes to nothing -- an empty sequence, or one of only gap
  characters -- is dropped when the file is read, rather than carried as far as the
  alignment and refused there. It produced no hits either way, but it used to print
  a pair header first:

  ```
  query 1: aso1 (20 nts) vs. target 1: empty (0 nts)
  ```

  Those headers are gone. Only the default and `-p1` output modes print them, so
  `-p2` and `-p3` are unaffected. Record numbers count every record in the file,
  including the dropped ones, so the surviving records keep the numbers they were
  reported under before.
