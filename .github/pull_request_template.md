## What changed

<!-- One paragraph. What does this PR make the app do that it did not do before? -->

## Scientific impact

<!-- REQUIRED for any change to the engine.
     What changes about the SEQUENCES the app produces? If nothing, say "none".
     If you added or changed a rule: what is the evidence, and what does enforcing
     it cost on the other objectives? -->

## Lane

<!-- Which lane owns this? Confirm you did not edit another lane's directory. -->

## Checklist

- [ ] I stayed inside my lane's directory
- [ ] I did not add a dependency or touch `uv.lock`
- [ ] Any new rule carries citations, an evidence badge and `last_verified`
- [ ] Any new SOFT rule explains its default weight in `weight_provenance`
- [ ] I did not weaken, skip or `xfail` an existing test
- [ ] `ruff check`, `mypy` and `pytest` pass locally
- [ ] If I added a CI job, I added it to `required-checks.needs`

## Owner sign-off

<!-- Load-bearing paths (the oracle, core/, invariants, benchmarks, data/, .github/)
     require the matching approved:* label before merge. -->
