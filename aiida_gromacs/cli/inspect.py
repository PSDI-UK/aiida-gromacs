#!/usr/bin/env python
"""
Inspect AiiDA databases.
"""

# from aiida import manage, orm, profile_context
from aiida.storage.sqlite_zip.backend import SqliteZipBackend

archive_profile = SqliteZipBackend.create_profile('test.aiida')
print(archive_profile)