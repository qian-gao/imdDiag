### -------------------------------- Import test data ---------------------- ### 

data.test <- read_xlsx("data/IMD_3_test_data.xlsx")

### -------------------------------- Start run ----------------------------- ### 

data.start <- read_xlsx("data/IMD_1_starting_data.xlsx")
ref.start <-
  variable_selection(
    input.control = data.start[data.start$IMD == "Control", ],
    input.sample = data.start[data.start$IMD != "Control", ],
    reference = NULL,
    start.run = TRUE)

### -------------------------------- Train more 1 --------------------------- ### 

data.more.1 <- read_xlsx("data/IMD_2_01_more_data.xlsx")
ref.more.1 <-
  variable_selection(
    input.control = data.more.1[data.more.1$IMD == "Control", ],
    input.sample = data.more.1[data.more.1$IMD != "Control", ],
    reference = ref.start,
    start.run = FALSE)

### -------------------------------- Train more 2 --------------------------- ### 

data.more.2 <- read_xlsx("data/IMD_2_02_more_data.xlsx")
ref.more.2 <-
  variable_selection(
    input.control = data.more.2[data.more.2$IMD == "Control", ],
    input.sample = data.more.2[data.more.2$IMD != "Control", ],
    reference = ref.more.1,
    start.run = FALSE)

### -------------------------------- Disgnosis run ------------------------- ### 

reference <- ref.more.2
reference <- ref.start

diagnosis <-
  generate_diagnosis( 
    input.control = data.test[data.test$IMD == "Control", ],
    input.sample = data.test[data.test$IMD != "Control", ],
    reference = reference, 
    top.nr = 3)




