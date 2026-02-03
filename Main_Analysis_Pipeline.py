import streamlit as st
import sqlite3
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

### SQL Database Scheme query

csv_path = r'cell-count.csv'
db_path = r'Loblaw_Bio_cell_count.db'
if os.path.exists(db_path):
    os.remove(db_path)
    
st.set_page_config(page_title='Loblaw Bio: Immune Cell Population Analysis')
st.title('Loblaw Bio: Clinical Trial Dashboard')

### Part 1 ================================================================================================================================


def initialize_database(csv_path, db_path):
	conn = sqlite3.connect(db_path)
	cursor = conn.cursor()
	cursor.execute('PRAGMA foreign_keys = ON;')

	df_raw = pd.read_csv(csv_path)
	df_raw.to_sql('cell-count', conn, if_exists='replace', index=False)

	sql_script = """
	/* =========================================================================================================
	DATA CHECK
	Preliminary checks to ensure data robustness. Commented out in final product to keep execution clean.
	=========================================================================================================
	*/

	/*
	-- Ensure Primary Key has no duplicates
	SELECT sample, COUNT(*) 
	FROM "cell-count" 
	GROUP BY sample
	HAVING COUNT(*) > 1;

	-- Ensure no NULLs in data for critical columns
	SELECT COUNT(*) AS missing_ids
	FROM "cell-count"
	WHERE subject IS NULL OR subject = '';
	*/

	/* 
	=========================================================================================================
	PRODUCTION SCHEMA 
		Table cleaning, separation, and normalization
	=========================================================================================================
	*/

	-- Subject Table Processing
	-- Define new table
	DROP TABLE IF EXISTS subject_table;
	CREATE TABLE subject_table (
		subject TEXT PRIMARY KEY,
		condition TEXT,
		age INTEGER,
		sex TEXT,
		treatment TEXT,
		response TEXT
	);

	-- Data insertion, select only distinct values, and normalize subject primary key via padding
	INSERT INTO subject_table (subject, condition, age, sex, treatment, response)
	SELECT DISTINCT 'sbj' || PRINTF('%04d', CAST(SUBSTR(subject, 4) AS INTEGER)), condition, age, sex, treatment, response
	FROM "cell-count";

	/*
	-- sanity check for duplicates
	SELECT 
		COUNT(*) AS total_rows, 
		COUNT(DISTINCT subject) AS unique_subjects
	FROM subject_table;
	*/


	-- Project Table Processing

	-- Check to see if project corresponds to one label each
	/*
	SELECT 
		project, 
		COUNT(DISTINCT project) AS label_count,
		GROUP_CONCAT(DISTINCT sample_type) AS labels_found
	FROM "cell-count"
	GROUP BY project;
	*/

	DROP TABLE IF EXISTS project_table;
	CREATE TABLE project_table (
		project TEXT PRIMARY KEY,
		sample_type TEXT
	);

	INSERT INTO project_table (project, sample_type)
	SELECT DISTINCT project, sample_type
	FROM "cell-count";


	-- Sample Table Processing

	-- Create and define sample table
	DROP TABLE IF EXISTS sample_table;
	CREATE TABLE sample_table (
		sample TEXT PRIMARY KEY,
		subject TEXT,
		project TEXT,
		time_from_treatment_start INTEGER,
		b_cell REAL,
		cd8_t_cell REAL,
		cd4_t_cell REAL,
		nk_cell REAL,
		monocyte REAL,
		FOREIGN KEY (subject) REFERENCES subject_table(subject),
		FOREIGN KEY (project) REFERENCES project_table(project)
	);

	INSERT INTO sample_table (sample, subject, project, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
	SELECT sample, 'sbj' || PRINTF('%04d', CAST(SUBSTR(subject, 4) AS INTEGER)), project, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte
	FROM "cell-count";
	"""

	cursor.executescript(sql_script)
	conn.commit()
	print(f'Database scheme successful')
	return conn

conn = initialize_database(csv_path, db_path)