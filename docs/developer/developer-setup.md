# Developer Setup

## Install uv

For a complete installation overview, go to the [UV GitHub repository](https://github.com/astral-sh/uv#installation)

```bash
# On macOS and Linux.
curl -LsSf https://astral.sh/uv/install.sh | sh
```

To upgrade uv in-place to the latest version, use :

```bash
uv self update
```

### Working with uv

- [Running commands in projects](https://docs.astral.sh/uv/concepts/projects/run/)
- [Managing dependencies](https://docs.astral.sh/uv/concepts/projects/dependencies/)
- [Configuring projects](https://docs.astral.sh/uv/concepts/projects/config/)

## Install prek

For a complete installation overview, go to the [prek GitHub repository](https://github.com/j178/prek#installation)

```bash
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/j178/prek/releases/download/v0.5.9/prek-installer.sh | sh
```

This installs prek version `v0.5.9`.

To upgrade prek in-place to the latest version, use :

```bash
prek self update
```

### Working with prek

- [Quickstart](https://prek.j178.dev/quickstart/)
- [Running hooks locally](https://prek.j178.dev/running-hooks/)
- [Configuration](https://prek.j178.dev/configuration/)
