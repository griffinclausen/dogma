init:
				pip install pipenv
				pipenv install --dev

test:
				pipenv run py.test tests

flake8:
				pipenv run flake8 --ignore=E501,F401,E128,E402,E731,F821 dogma
