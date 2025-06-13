Here is the source code of paper *Skyline Keyword-Aware Community Search on Semantic and Structure*

### Example of compiling and running (Linux)

Before running, you should:
+ Prepare the dataset and query files according to the templates
+ Check file names and file paths

Then run the following commands,

```shell
$ make
$ main "bk" "0" "3" "l" "3"
# For running queries, five parameters should be specified:
# "bk" is the name of dataset
# "0" means that the similarity threshold
# "3" means that we are running queries on varying experiments
# "l" is the type of attribute, where "i" is the textual attribute, and "l" is the locational attribute 
# "3" is the k value
```


