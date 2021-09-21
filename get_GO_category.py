#!/usr/bin/env python3
# coding: utf-8

import sys
import requests


'''Searchs GO categories (BP, MF, CC) for GO IDs using QuickGO API'''


def getGOcategory(goId):

    """
    Search GO category (BP, MF, CC) for a GO ID using QuickGO API

    Parameter: GO ID (ex: GO:0015979)
    
    Return: GO category (BP, MF or CC)
    """

    goId = goId.replace(':', '%3A')

    # search GO ID in QuickGO 
    requestURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/' + goId
    r = requests.get(requestURL, headers={ 'Accept' : 'application/json'})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # find GO ID category
    responseBody = r.text
    category = responseBody.split(':')[-3].split('\"')[1]

    # convert QuickGO category to MF, BP and CC format
    if category == 'molecular_function':
        category = 'MF'
    elif category == 'biological_process':
        category = 'BP'
    else:
        category = 'CC'

    return category

getGOcategory(sys.argv[1])