#!/usr/bin/env python
import csv
import json
import os
import re
from argparse import ArgumentParser
from collections import defaultdict

import atexit
from unittest import TestCase

import pandas as pd
import psycopg2
import requests
from ebi_eva_common_pyutils.config_utils import get_pg_metadata_uri_for_eva_profile
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.mongo_utils import get_mongo_connection_handle
from ebi_eva_common_pyutils.pg_utils import execute_query, get_pg_connection_handle, get_all_results_for_query

from tasks.eva_2406.gather_release_species import retrieve_assembly_summary_from_species_name, best_assembly, \
    most_recent_assembly, format_assembly_summaries

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


class TestGatherReleaseSpecies(TestCase):

    def test_retrieve_assembly_summary_from_species_name(self):
        assembly_list = (retrieve_assembly_summary_from_species_name('Bos taurus'))
        format_assembly_summaries([most_recent_assembly(assembly_list)])
        format_assembly_summaries([best_assembly(assembly_list)])