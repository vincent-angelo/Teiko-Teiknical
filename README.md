# Teiko_Loblaw_Bio

## How to run
Dependencies are listed in requirements.txt. They can be installed by running ```pip install -r requirements.txt``` in the terminal.
Please run the code through ```streamlit run Analysis.py```

Dashboard link: https://teikoloblawbio-wbwtpe8un8du2xwf7usoy9.streamlit.app/

## Relational Database Schema and Rationale
Given the cell-count.csv file, the raw data is separated into three tables through SQL: 
- subject_table: Stores unique patient metadata (subject [Primary Key], Condition, Age, Sex, Treatment, Response).
- project_table: Stores unique project identifiers and associated sample types (project [Primary Key], sample type).
- sample_table: Stores individual sample observations, including timepoints and cell counts. It uses Foreign Keys to link back to the parent tables.

This scheme has normalized tables that fulfill the requirements of the third normal form (3NF) tables. As a result, the scheme avoids redundancy, ensuring that any data is stored only once and is connected between tables through foreign keys, linking every sample record to its respective subject and project.

This efficient layout scales well, as 3NF structures help minimize data storage requirements by eliminating the repetition of metadata, leaving more space for high-volume sample observations. Furthermore, this allows any data correction to happen at one point, instead of multiple, such as a correction on subject data would only require an edit in the subject_table once, instead of multiple times in the original table. Lastly, should the company introduce additional products/services, the new data can simply be tied to the subject data to keep consistent data between services. 

### Notes on data management rationale:
- 'sbj' was purposefully left behind and not cleaned (sbj0001 -> 0001) with scaling in mind. To account for more data entries in the future, where 'sbj' is not guaranteed and clinicians can assign a different prefix to the subjects, 'sbj' was not removed to avoid ambiguity and allow such changes when scaling up.
- Subject values are padded (sbj001 -> sbj0001) to ensure alphanumeric sorting and create a more intuitive table organization for clinicians/clients to review. Depending on the scale of the data, the amount of padding could be easily or dynamically adjusted to accommodate the larger number of subjects.
- Sanity and data integrity checks are provided but commented out in the query to ensure efficient execution.
- Although the current analysis retrieves data from all tables for comprehensive analysis, when the number of tables and samples increases in the future, this scheme is still optimal to enable selective retrieval and reduces memory overhead compared to processing the raw table. 

Given the flexibility of the modified database schema, there are a few interesting analyses that might be possible as the data sama sample size scales up, such as:
- Changes in cell count over time as new treatment is introduced at various points in time.
- Cross-comparison between old project data and newer projects covering different treatments or services.
- Training a predictive model to predict response based on cell counts and other newly gathered features.

## Code Structure

The code features a hybrid approach, where the relational database is built through SQL prompts, while the table retrieval and data analysis are done through Python (pandas), and the interface design was done through Streamlit. 

- Data Management (SQL): The cell-count.csv raw data was processed by a SQL prompt, separating it into three tables (mentioned above) and storing them in a .db file to help enforce Relational Integrity through Foreign Key constraints. SQL provides an efficient storage and retrieval process, which will be crucial as the database size scales up. 
- Analytical Layer (Python/Pandas): Afterward, for each part of the problem given, I retrieved the necessary and relevant data from the database and converted each to a pandas dataframe. Further processing is done through pandas in Python to tailor to the needs of every question. The code is dynamically coded to ensure that further additions or edits would be a simple change, such as by having my analysis logic (melting, grouping, and calculating frequencies) driven by the column headers and unique values present in the database. The boxplot for part 3 is created through seaborn, and the significance test is done through Welch's T-test that iterates through all available cell populations, setting the significance threshold at <0.05, a commonly used value. 
- Presentation (Streamlit): Lastly, all the analysis was wrapped into a presentable format for Yah D’yada through Streamlit. This helps create an interactive dashboard to organize and show all the analysis results, highlighting only the important results without all the analysis clutter. 

