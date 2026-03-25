# Container Build

This project expects a Boltz runtime image at `containers/current`.

On `cbe`, the recommended production build is an Apptainer sandbox created from
the existing Boltz image with Plotly and Kaleido baked in:

```bash
apptainer build --fakeroot --sandbox containers/boltz_screen containers/boltz_screen.def
ln -sfn boltz_screen containers/current
```

If your Apptainer installation supports it reliably, you can convert the tested
sandbox into an immutable SIF later and repoint `containers/current` to that
file instead.
