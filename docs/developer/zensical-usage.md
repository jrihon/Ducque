# Zensical Usage

The documentation is rendered with [Zensical](https://zensical.org/).
It is installed as part of the `dev` dependency group, so no separate installation is needed after `uv sync`.

The configuration lives in `zensical.toml` at the root of the repository.

## Serve the documentation locally

Run the following from the root of the repository:

```bash
uv run zensical serve
```

This builds the documentation and serves it on [localhost:8000](http://localhost:8000), rebuilding whenever a file changes.
