
# Template repository for a library project

## Initial set up

```bash
pixi run init
```

## Incomporate updates changes to upstream fork

If you forked the `munch-group-project` rather than using it as template, you can incorporate changes/fixes made to `munch-group-project`.

Add upstream if not already added

```bash
git remote add upstream https://github.com/munch-group/munch-group-project.git
```

Fetch upstream changes

```bash
git fetch upstream
```

Either rebase your changes on top of upstream (cleaner history)

```bash
git rebase upstream/main
```

Or, merge upstream into your fork (preserves history)

```bash
git merge upstream/main
```

If you want to see what's changed upstream before applying:

```bash
git log HEAD..upstream/main
```

See the actual diff

```bash
git diff HEAD...upstream/main
```

Then push your updated fork:

```bash
git push origin main
```

If you rebased and need to force push
    
```bash
git push origin main --force-with-lease
```
