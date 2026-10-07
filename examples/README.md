# Examples

Sample input files for manual testing and for CI.

`params.json` holds the `--param`/`--params-json` values `cauldron job run` uses against this
plugin's inputs (see `plugin.yaml`'s `inputs:` section, or run
`cauldron plugin inputs lysoip-gating-ranking` once installed). CI runs this automatically via
`.github/workflows/test-plugin.yml`.

```json
{
  "differential_expression_file": "examples/differential_expression.tsv",
  "reproducibility_file": "examples/reproducibility.tsv",
  "organelle_specificity_file": "examples/organelle_specificity.tsv",
  "go_evidence_file": "examples/go_evidence.tsv",
  "compass_file": "examples/compass.tsv",
  "crapome_file": "examples/crapome.tsv"
}
```
