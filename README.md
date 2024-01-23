# Simple Cartesian Mesh

# Build

Currently, this project depends on an unreleased version of _flecsi_.
You can build using the following spack environment:

```python
spack:
  specs:
  - flecsi@2.3-beta
  - yaml-cpp
  view: true
  concretizer:
    unify: true
  repos:
  - /PATH_TO_PROJECT/fdm/default/spack-repo
```
where _PATH\_TO\_PROJECT_ is set to your git clone.

<!-- vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 : -->
