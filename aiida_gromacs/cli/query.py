#!/usr/bin/env python
"""
Query AiiDA databases.
"""

from aiida.orm.querybuilder import QueryBuilder
from aiida.plugins import CalculationFactory
from aiida import load_profile
from aiida.manage.configuration import get_profile

profile = load_profile()


MyAppCalculation = CalculationFactory("general-MD")

qb = QueryBuilder()
qb.append(MyAppCalculation) #, tag="calcjob")
first_row = qb.first()
all_results_d = qb.dict()
print(all_results_d)

node_uuids = [node[0].uuid for node in qb.all()]
print(node_uuids)

first_row[0].get_incoming()




# archive = Archive()
# archive.add_nodes(node_uuids)
# archive.archive('aiida-data.zip')