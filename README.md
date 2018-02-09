# Templates

A bunch of templates in python for processing HTS data. Particularly
useful for using with nextflow pipelines.

## Quick reference

* mapping2json.py - exports results from a samtool depth file to a json
file that contains a `key:value` such as `accession number:coverage` .

* mashdist2json.py - exports results from `mash dist` to a json file
that contains a `key:value` such as `accession number:distance` .

* mashscreen2json.py - exports results from `mash screen` to a json
file that contains a `key:[values]` such as `accession number:[copy number, identity]` .

## How to use as a submodule

### Add templates to your project

```
git submodule add https://github.com/ODiogoSilva/templates.git templates
```

### Update templates on your project

```
git submodule foreach git pull origin master
```