PYTHON ?= python
OUT    ?= outputs/figures

.PHONY: all install validate figures checksums clean

all: validate figures

install:
	$(PYTHON) -m pip install -r requirements.txt

validate:
	$(PYTHON) scripts/validate_bundle.py --report outputs/validation.json

figures:
	$(PYTHON) scripts/reproduce_figures.py --out $(OUT)

checksums:
	$(PYTHON) scripts/update_checksums.py

clean:
	rm -rf outputs
