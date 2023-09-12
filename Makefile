MODULE := simulations
BLUE='\033[0;34m'
NC='\033[0m' 
TAG := 1 #$(shell git describe --tags --always --dirty)f

run_multiprocessing:
ifneq ($(and $(conf),$(n_models)),)
	@python3 -m $(MODULE) 0 $${conf} $${n_models}
else 
	@echo "Usage: make run_multiprocessing conf=conf_name.yml n_models=10 \nThis simulates 10 models as specified in conf_name.yml \nConfiguration file should be save under in /config/conf_name.yml"
endif

run_optuna:
ifneq ($(and $(setting), $(name),$(conf),$(n_trials)),)
	@python3 -m $(MODULE) $${setting} $${name} $${conf} $${n_trials}
else 
	@echo "Usage: make run_optuna setting=n name=name conf=conf_name.yml n_trials=100 \n"
endif

test:
	@pytest tests

lint:
	@echo "\n${BLUE}Running Pylint against source and test files...${NC}\n"
	@pylint --rcfile=setup.cfg **/*.py
	@echo "\n${BLUE}Running Flake8 against source and test files...${NC}\n"
	@flake8

clean:
	rm -rf .pytest_cache .coverage .pytest_cache coverage.xml tmp

version:
	@echo $(TAG)

.PHONY: clean test