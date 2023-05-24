#!/usr/bin/env python
"""
Query AiiDA databases.
"""

from aiida.orm import QueryBuilder


MyAppCalculation = CalculationFactory("general-MD")

qb = QueryBuilder()
qb.append(MyAppCalculation, tag="calcjob")