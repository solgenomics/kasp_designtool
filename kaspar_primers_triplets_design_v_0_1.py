#! /usr/bin/python -O
# -*- coding: utf-8 -*-
# kaspar_primers_triplets_design_v_0_1.py
__version__ = "0.1"

"""kaspar_primers_triplets_design_v_0_1.py
This program designs one or several triplets of primers from a list of SNPs
"""

# imported modules
# import commands, os, string, sys, csv, shutil
import commands, os, string, sys, csv, shutil

# list of functions
# chop(s)
# error(no, s = "")
# files_close()
# help()
# help_detailed()
# getting_default_pcr_parameters()
# primers_triplets_design()
# writing_primer3_input_files()
# testing_sequence(seq,seq_id,snp_nb)
# primer3_parser(fl_primer3_output_name,fl_primer3_output_parser_name,direction)
# primers_triplets_file_writing(fl_fam_primer3_output_parser_name,fl_hex_primer3_output_parser_name)
# writing_all_primers_triplets_files()
# testing_primers_triplets_specificities()
# primers_triplet_specificity(i,j,direction)
# adding_specificities_informations()
# writing_primers_triplets_specificities_file()
# print_msg(msg = "\n")



#======================================================#
# function chop(s) : trimming \n at the end of a string #
#======================================================#

def chop(s) :

  if type(s) == type('a') :

    if s[-1] == '\n' :
      s = s[:-1]

  return s



#======================================#
# function error() : to process errors #
#======================================#

def error(no, s = "") :

  lst_errors = [ \
    (1, "Non integer argument used in error call (", ").", 1),
    (2, "Unknown error number (", ").", 1),
    (3, "Wrong list of arguments ", "(use -h or -H argument to get information).", 1),
    (4, "Mandatory argument ", " is missing (use -h or -H argument to get information).", 1),
    (5, "Impossible to open file '", "' with mode 'a (append)'.", 1),
    (6, "Impossible to open file '", "' with mode 'r (read)'.", 1),
    (7, "Impossible to open file '", "' with mode 'w (write)'.", 1),
    (8, "The directory '", "' doesn't exist.", 1)]

  no_error = 0
  label   = ""

  if type(no) == type(1) :

    if no < 1 or no > len(lst_errors) :
      no_error = 2
      label   = str(no)

    else :
      no_error = no
      label   = s

  else :
    no_error = 1
    label   = no

  t_error = lst_errors[no_error - 1]

  if fl_err :
    fl_err.write("\nError %d :\n" % (no_error))
    fl_err.write("%s\n" % (t_error[1] + label + t_error[2]))
    fl_err.flush()
  else :
    sys.stderr.write("\nError %d :\n" % (no_error))
    sys.stderr.write("%s\n" % (t_error[1] + label + t_error[2]))
    sys.stderr.flush()

  if t_error[3] == 0 :
    return

  sys.exit(1)



#====================================================#
# function files_close() : to close all opened files #
#====================================================#

def files_close() :

  if fl_input :
    fl_input.close()

  if fl_pcr_parameters :
    fl_pcr_parameters.close()

  if fl_err :
    fl_err.close ()

  return



#============================================================#
# function help() : to display information about the program #
#============================================================#

def help() :

  print """kaspar_primers_triplets_design_v_0_1.py [-h] [-H] [-pcr_parameters_default] -in f1 [-flanking_sequence_min_size l1] [-checked_zone_length l2] [-pcr_parameters_file f2] [-pcr_parameters p1 p2 ...] [-species_db db1] [-err f3] [-debug]
  This program designs several triplets of primers.
  Arguments:
  -h                             : to get some information about this program (what you get now)
  -H                             : to get detailed information about this program
  -pcr_parameters_default        : to get the PCR parameters by default
  -in f1                         : file containing a list of SNP with their backward and forward sequences
  -flanking_sequence_min_size l1 : minimum length of the flanking sequences
  -checked_zone_length l2        : no IUPAC code among the l2 nucleotides before and after the SNP in a sequence
  -pcr_parameters_file f2        : file containing the PCR parameters used instead of the default PCR parameters
  -pcr_parameters p1 p2 ...      : PCR parameters used instead of the default PCR parameters
  -species_db db1                : path to the data base of the species used to check the specificities of the triplets of primers
  -err f3                        : errors file (if not defined standard error is used)
  -debug                         : debug mode
 """
  sys.exit(0)



#=====================================================================#
# function help_detailed() : to display information about the program #
#=====================================================================#

def help_detailed() :

  help()



#======================================================================================#
# function getting_default_pcr_parameters() : to display the PCR parameters by default #
#======================================================================================#

def getting_default_pcr_parameters() :

  line = ""

  # fl_pcr_parameters_default_name : file containing the default PCR parameters
  try :
    fl_pcr_parameters_default = open(fl_pcr_parameters_default_name, 'r')
  except :
    error(6, fl_pcr_parameters_default_name)

  while 1 :
    line = fl_pcr_parameters_default.readline()

    if line == "" :
      # EOF ...
      break

    print line[:-1]

  if fl_pcr_parameters_default :
    fl_pcr_parameters_default.close()

  sys.exit(0)



#=================================================================================================#
# function primers_triplets_design() : to design one or several triplets of primers for each SNP  #
#=================================================================================================#

def primers_triplets_design() :

  line      = ""
  lst_field = []

  # creation of the directory SNP in the current working directory
  if os.path.exists(current_working_directory + "SNP") :
    shutil.rmtree(current_working_directory + "SNP")
  os.mkdir(current_working_directory + "SNP")

  # writing the Primer3 input files
  writing_primer3_input_files()

  os.chdir(current_working_directory + "SNP")
  # file including the number of triplets by SNP
  fl_nb_of_triplets_by_snp = open("nb_of_triplets_by_snp.csv",'w')
  fl_nb_of_triplets_by_snp.write("id;total;forward;reverse \n")

  i = 1
  while i <= nb_of_snp :

    # the SNP is processed
    if snp_processing[i-1] == 1 :

      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/FORWARD")
      # executing Primer3 for the FAM and HEX alleles of the SNP
      os.system("/usr/bin/primer3_core -output=fam_primer3_output.txt < fam_primer3_input.txt")
      os.system("/usr/bin/primer3_core -output=hex_primer3_output.txt < hex_primer3_input.txt")
      # parsing the Primer3 output files
      primer3_parser("fam_primer3_output.txt","fam_primer3_output_parsed.csv","FORWARD")
      primer3_parser("hex_primer3_output.txt","hex_primer3_output_parsed.csv","FORWARD")
      #writing the file containing the triplets of primers from the 2 previous files
      primers_triplets_file_writing("fam_primer3_output_parsed.csv","hex_primer3_output_parsed.csv")

      forward_nb_of_triplets.append(nb_of_triplets)  # nb_of_triplets is calculated in the primers_triplets_file_writing function

      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/REVERSE")
      # executing Primer3 for the FAM and HEX alleles of the SNP
      os.system("/usr/bin/primer3_core -output=fam_primer3_output.txt < fam_primer3_input.txt")
      os.system("/usr/bin/primer3_core -output=hex_primer3_output.txt < hex_primer3_input.txt")
      # parsing the Primer3 output files
      primer3_parser("fam_primer3_output.txt","fam_primer3_output_parsed.csv","REVERSE")
      primer3_parser("hex_primer3_output.txt","hex_primer3_output_parsed.csv","REVERSE")
      # writing the file containing the triplets of primers from the 2 previous files
      primers_triplets_file_writing("fam_primer3_output_parsed.csv","hex_primer3_output_parsed.csv")

      reverse_nb_of_triplets.append(nb_of_triplets)  # nb_of_triplets is calculated in the primers_triplets_file_writing function

      os.chdir(current_working_directory + "SNP")
      total_nb_of_triplets = forward_nb_of_triplets[i-1] + reverse_nb_of_triplets[i-1]
      fl_nb_of_triplets_by_snp.write(lst_snp_id[i-1] + ";" + str(total_nb_of_triplets) + ";" + str(forward_nb_of_triplets[i-1]) + ";" + str(reverse_nb_of_triplets[i-1]) + "\n")

    else :

      forward_nb_of_triplets.append(0)
      reverse_nb_of_triplets.append(0)
      fl_nb_of_triplets_by_snp.write(lst_snp_id[i-1] + ";not processed \n")

    i += 1

  if fl_nb_of_triplets_by_snp :
    fl_nb_of_triplets_by_snp.close()

  # testing the specificities of the triplets of primers
  if len(species_db_name) :
    testing_primers_triplets_specificities()

  # writing files including all the triplets of primers
  writing_all_primers_triplets_files()



#===========================================================================#
# function writing_primer3_input_files() : to write the Primer3 input files #
#===========================================================================#

def writing_primer3_input_files() :

  global nb_of_snp
  line      = ""
  lst_field = []

  while 1 :
    line = fl_input.readline()

    if line == "" :
        # EOF ...
        break

    if line.find(">") == 0 :

      nb_of_snp = nb_of_snp + 1

      lst_field = line.split('>')
      snp_id = lst_field[1]
      lst_snp_id.append(chop(snp_id))

      try:
        os.mkdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp))
      except OSError:
        pass

      try:
        os.mkdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
      except OSError:
        pass

      try:
        os.mkdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
      except OSError:
        pass

      line = fl_input.readline()

      # testing the sequence
      testing_sequence(line,snp_id,nb_of_snp)

       # the SNP is processed
      if snp_processing[nb_of_snp-1] == 1 :

        lst_field = line.split('[')
        lst_field_02 = lst_field[1].split(']')

        snp_position = len(lst_field[0])

        if lst_field_02[0].find('/') != -1 :
          lst_field_03 = lst_field_02[0].split("/")
          sequence_template_01 = lst_field[0] + lst_field_03[0] + lst_field_02[1]
          sequence_template_02 = lst_field[0] + lst_field_03[-1] + lst_field_02[1]
        elif lst_field_02[0] == 'R' :
          sequence_template_01 = lst_field[0] + 'A' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'G' + lst_field_02[1]
        elif lst_field_02[0] == 'Y' :
          sequence_template_01 = lst_field[0] + 'C' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'T' + lst_field_02[1]
        elif lst_field_02[0] == 'S' :
          sequence_template_01 = lst_field[0] + 'G' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'C' + lst_field_02[1]
        elif lst_field_02[0] == 'W' :
          sequence_template_01 = lst_field[0] + 'A' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'T' + lst_field_02[1]
        elif lst_field_02[0] == 'K' :
          sequence_template_01 = lst_field[0] + 'G' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'T' + lst_field_02[1]
        elif lst_field_02[0] == 'M' :
          sequence_template_01 = lst_field[0] + 'A' + lst_field_02[1]
          sequence_template_02 = lst_field[0] + 'C' + lst_field_02[1]

        os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
        fl_allele_fam = open("fam_primer3_input.txt", 'w')
        fl_allele_hex = open("hex_primer3_input.txt", 'w')
        fl_allele_fam.write("SEQUENCE_ID=" + snp_id)
        fl_allele_fam.write("SEQUENCE_TEMPLATE=" + sequence_template_01)
        fl_allele_fam.write("SEQUENCE_FORCE_LEFT_END=" + str(snp_position) + "\n")
        fl_allele_hex.write("SEQUENCE_ID=" + snp_id)
        fl_allele_hex.write("SEQUENCE_TEMPLATE=" + sequence_template_02)
        fl_allele_hex.write("SEQUENCE_FORCE_LEFT_END=" + str(snp_position) + "\n")
        if fl_allele_fam :
          fl_allele_fam.close()
        if fl_allele_hex :
          fl_allele_hex.close()

        os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
        fl_allele_fam = open("fam_primer3_input.txt", 'w')
        fl_allele_hex = open("hex_primer3_input.txt", 'w')
        fl_allele_fam.write("SEQUENCE_ID=" + snp_id)
        fl_allele_fam.write("SEQUENCE_TEMPLATE=" + sequence_template_01)
        fl_allele_fam.write("SEQUENCE_FORCE_RIGHT_END=" + str(snp_position) + "\n")
        fl_allele_hex.write("SEQUENCE_ID=" + snp_id)
        fl_allele_hex.write("SEQUENCE_TEMPLATE=" + sequence_template_02)
        fl_allele_hex.write("SEQUENCE_FORCE_RIGHT_END=" + str(snp_position) + "\n")
        if fl_allele_fam :
          fl_allele_fam.close()
        if fl_allele_hex :
          fl_allele_hex.close()

        # fl_pcr_parameters_name : file containing the PCR parameters
        if len(fl_pcr_parameters_name) :

          try :
            fl_pcr_parameters = open(fl_pcr_parameters_name, 'r')
          except :
            error(6, fl_pcr_parameters_name)

          line_pcr_parameters = ""

          while 1 :
            line_pcr_parameters = fl_pcr_parameters.readline()

            if line_pcr_parameters == "" :
              # EOF ...
              break

            os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
            fl_allele_fam = open("fam_primer3_input.txt", 'a')
            fl_allele_hex = open("hex_primer3_input.txt", 'a')
            fl_allele_fam.write(line_pcr_parameters)
            fl_allele_hex.write(line_pcr_parameters)
            if fl_allele_fam :
              fl_allele_fam.close()
            if fl_allele_hex :
              fl_allele_hex.close()

            os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
            fl_allele_fam = open("fam_primer3_input.txt", 'a')
            fl_allele_hex = open("hex_primer3_input.txt", 'a')
            fl_allele_fam.write(line_pcr_parameters)
            fl_allele_hex.write(line_pcr_parameters)
            if fl_allele_fam :
              fl_allele_fam.close()
            if fl_allele_hex :
              fl_allele_hex.close()

          if fl_pcr_parameters :
            fl_pcr_parameters.close()

        else :

          # fl_pcr_parameters_default_name : file containing the default PCR parameters
          try :
            fl_pcr_parameters_default = open(fl_pcr_parameters_default_name, 'r')
          except :
            error(6, fl_pcr_parameters_default_name)

          line_pcr_parameters_default = ""

          user_modif_param = 0

          while 1 :
            line_pcr_parameters_default = fl_pcr_parameters_default.readline()

            if line_pcr_parameters_default == "" :
              # EOF ...
              break

            j = 0

            while j < len(list_pcr_parameters) :

              line_pcr_parameters_default_split = line_pcr_parameters_default.split("=")
              list_pcr_parameters_split = list_pcr_parameters[j].split("=")

              user_modif_param = 0

              # the default PCR parameter is modified by the user (new value given in argument)
              if line_pcr_parameters_default_split[0] == list_pcr_parameters_split[0] :

                os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
                fl_allele_fam = open("fam_primer3_input.txt", 'a')
                fl_allele_hex = open("hex_primer3_input.txt", 'a')
                fl_allele_fam.write(list_pcr_parameters[j] + "\n")
                fl_allele_hex.write(list_pcr_parameters[j] + "\n")
                if fl_allele_fam :
                  fl_allele_fam.close()
                if fl_allele_hex :
                  fl_allele_hex.close()
                os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
                fl_allele_fam = open("fam_primer3_input.txt", 'a')
                fl_allele_hex = open("hex_primer3_input.txt", 'a')
                fl_allele_fam.write(list_pcr_parameters[j] + "\n")
                fl_allele_hex.write(list_pcr_parameters[j] + "\n")
                if fl_allele_fam :
                  fl_allele_fam.close()
                if fl_allele_hex :
                  fl_allele_hex.close()

                user_modif_param = 1
                break

              j += 1

            # the default PCR parameter is not modified by the user; the default value is written in the FAM and HEX alleles files
            if user_modif_param == 0 :
              os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
              fl_allele_fam = open("fam_primer3_input.txt", 'a')
              fl_allele_hex = open("hex_primer3_input.txt", 'a')
              fl_allele_fam.write(line_pcr_parameters_default)
              fl_allele_hex.write(line_pcr_parameters_default)
              if fl_allele_fam :
                fl_allele_fam.close()
              if fl_allele_hex :
                fl_allele_hex.close()
              os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
              fl_allele_fam = open("fam_primer3_input.txt", 'a')
              fl_allele_hex = open("hex_primer3_input.txt", 'a')
              fl_allele_fam.write(line_pcr_parameters_default)
              fl_allele_hex.write(line_pcr_parameters_default)
              if fl_allele_fam :
                fl_allele_fam.close()
              if fl_allele_hex :
                fl_allele_hex.close()

          if fl_pcr_parameters_default :
            fl_pcr_parameters_default.close()

        # adding "=" at the end of allele_fam.txt and allele_txt
        os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/FORWARD")
        fl_allele_fam = open("fam_primer3_input.txt", 'a')
        fl_allele_hex = open("hex_primer3_input.txt", 'a')
        fl_allele_fam.write("=")
        fl_allele_hex.write("=")
        if fl_allele_fam :
          fl_allele_fam.close()
        if fl_allele_hex :
          fl_allele_hex.close()
        os.chdir(current_working_directory + "SNP/SNP_" + str(nb_of_snp) + "/REVERSE")
        fl_allele_fam = open("fam_primer3_input.txt", 'a')
        fl_allele_hex = open("hex_primer3_input.txt", 'a')
        fl_allele_fam.write("=")
        fl_allele_hex.write("=")
        if fl_allele_fam :
          fl_allele_fam.close()
        if fl_allele_hex :
          fl_allele_hex.close()



#==============================================================================#
# function testing_sequence(seq,seq_id,snp_nb) : to test the sequence of a SNP #
#==============================================================================#

def testing_sequence(seq,seq_id,snp_nb) :

  # no SNP is indicated in the sequence
  if seq.count('[') == 0 and seq.count(']') == 0 :

    message = "in sequence " + str(seq_id[:-1]) + ", no SNP is indicated; this sequence wont' be processed"
    print_msg(message)

    snp_processing.append(0) # the sequence won't be processed

  # more than one SNP are indicated in the sequence
  elif seq.count('[') > 1 and seq.count('[') == seq.count(']') :

    message = "in sequence " + str(seq_id[:-1]) + ", more than one SNP are indicated; this sequence won't be processed"
    print_msg(message)

    snp_processing.append(0) # the sequence won't be processed

  # one (and only one) SNP is indicated in the sequence
  else :

    left_fl_seq = seq.split('[')[0]
    right_fl_seq = seq.split('[')[1].split(']')[1]

    # testing the flanking sequences
    if len(left_fl_seq) >= flanking_sequence_min_size and len(right_fl_seq)-1 >= flanking_sequence_min_size and \
       left_fl_seq[-checked_zone_length:].count('R') == 0 and left_fl_seq[-checked_zone_length:].count('Y') == 0 and left_fl_seq[-checked_zone_length:].count('S') == 0 and \
       left_fl_seq[-checked_zone_length:].count('W') == 0 and left_fl_seq[-checked_zone_length:].count('K') == 0 and left_fl_seq[-checked_zone_length:].count('M') == 0 and \
       right_fl_seq[:checked_zone_length].count('R') == 0 and right_fl_seq[:checked_zone_length].count('Y') == 0 and right_fl_seq[:checked_zone_length].count('S') == 0 and \
       right_fl_seq[:checked_zone_length].count('W') == 0 and right_fl_seq[:checked_zone_length].count('K') == 0 and right_fl_seq[:checked_zone_length].count('M') == 0 and \
       seq.count('[') == 1 and seq.count(']') == 1 :

      snp_processing.append(1) # the SNP will be processed

    else :

      snp_processing.append(0) # the SNP won't be processed

      if len(left_fl_seq) < flanking_sequence_min_size :
        message = "in sequence " + str(seq_id[:-1]) + ", left flanking sequence's length is lower than " + str(flanking_sequence_min_size) + " : SNP won't be processed"
        print_msg(message)

      if len(right_fl_seq)-1 < flanking_sequence_min_size :
        message = "in sequence " + str(seq_id[:-1]) + ", right flanking sequence's length is lower than " + str(flanking_sequence_min_size) + " : SNP won't be processed"
        print_msg(message)

      if left_fl_seq[-checked_zone_length:].count('R') > 0 or left_fl_seq[-checked_zone_length:].count('Y') > 0 or left_fl_seq[-checked_zone_length:].count('S') > 0 or \
         left_fl_seq[-checked_zone_length:].count('W') > 0 or left_fl_seq[-checked_zone_length:].count('K') > 0 or left_fl_seq[-checked_zone_length:].count('M') > 0 :

        message = "in sequence " + str(seq_id[:-1]) + ", there are one or more IUPAC code among the " + str(checked_zone_length) + " nucleotides before the SNP : SNP won't be processed"
        print_msg(message)

      if right_fl_seq[:checked_zone_length].count('R') > 0 or right_fl_seq[:checked_zone_length].count('Y') > 0 or right_fl_seq[:checked_zone_length].count('S') > 0 or \
         right_fl_seq[:checked_zone_length].count('W') > 0 or right_fl_seq[:checked_zone_length].count('K') > 0 or right_fl_seq[:checked_zone_length].count('M') > 0 :

        message = "in sequence " + str(seq_id[:-1]) + ", there are one or more IUPAC code among the " + str(checked_zone_length) + " nucleotides after the SNP : SNP won't be processed"
        print_msg(message)



#==========================================================================================================================#
# function primer3_parser(fl_primer3_output_name,fl_primer3_output_parser_name,direction) : to parse a Primer3 output file #
#==========================================================================================================================#

def primer3_parser(fl_primer3_output_name,fl_primer3_output_parser_name,direction) :

  fl_input = open(fl_primer3_output_name, 'r')
  fl_output = open(fl_primer3_output_parser_name, 'w')

  line      = ""
  lst_field = []

  fl_output.write("id;oligo;common primer;start;len;tm;%gc;any;3prime;specific primer;start;len;tm;%gc;any;3prime;pair's score \n")

  oligo_number = 0

  while 1 :
    line = fl_input.readline()

    if line == "" :
      # EOF ...
      break

    if line.find("SEQUENCE_ID") == 0 :
      sequence_id = chop(line.split("=")[1])

    if line.find("PRIMER_PAIR_" + str(oligo_number) + "_PENALTY") == 0 :
      pair_penalty = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "_SEQUENCE") == 0 :
      if direction == "FORWARD" :
        specific_sequence = chop(line.split("=")[1])
      else :
        common_sequence = chop(line.split("=")[1])

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "_SEQUENCE") == 0 :
      if direction == "FORWARD" :
        common_sequence = chop(line.split("=")[1])
      else :
        specific_sequence = chop(line.split("=")[1])

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "=") == 0 :
      lst_field_01 = line.split("=")
      lst_field_02 = lst_field_01[1].split(",")
      if direction == "FORWARD" :
        specific_start = chop(lst_field_02[0])
        specific_length = chop(lst_field_02[1])
      else :
        common_start = chop(lst_field_02[0])
        common_length = chop(lst_field_02[1])

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "=") == 0 :
      lst_field_01 = line.split("=")
      lst_field_02 = lst_field_01[1].split(",")
      if direction == "FORWARD" :
        common_start = chop(lst_field_02[0])
        common_length = chop(lst_field_02[1])
      else :
        specific_start = chop(lst_field_02[0])
        specific_length = chop(lst_field_02[1])

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "_TM") == 0 :
      if direction == "FORWARD" :
        specific_tm = round(float(chop(line.split("=")[1])),3)
      else :
        common_tm = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "_TM") == 0 :
      if direction == "FORWARD" :
        common_tm = round(float(chop(line.split("=")[1])),3)
      else :
        specific_tm = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "_GC_PERCENT") == 0 :
      if direction == "FORWARD" :
        specific_gc = round(float(chop(line.split("=")[1])),3)
      else :
        common_gc = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "_GC_PERCENT") == 0 :
      if direction == "FORWARD" :
        common_gc = round(float(chop(line.split("=")[1])),3)
      else :
        specific_gc = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "_SELF_ANY") == 0 :
      if direction == "FORWARD" :
        specific_any = round(float(chop(line.split("=")[1])),3)
      else :
        common_any = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "_SELF_ANY") == 0 :
      if direction == "FORWARD" :
        common_any = round(float(chop(line.split("=")[1])),3)
      else :
        specific_any = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_LEFT_" + str(oligo_number) + "_SELF_END") == 0 :
      if direction == "FORWARD" :
        specific_end = round(float(chop(line.split("=")[1])),3)
      else :
        common_end = round(float(chop(line.split("=")[1])),3)

    if line.find("PRIMER_RIGHT_" + str(oligo_number) + "_SELF_END") == 0 :
      if direction == "FORWARD" :
        common_end = round(float(chop(line.split("=")[1])),3)
      else :
        specific_end = round(float(chop(line.split("=")[1])),3)

      fl_output.write(sequence_id + ";OLIGO_" + str(oligo_number) + ";" + \
                        common_sequence + ";" + common_start + ";" + common_length + ";" + str(common_tm) + ";" + str(common_gc) + ";" + str(common_any) + ";" + str(common_end) + ";" + \
                        specific_sequence + ";" + specific_start + ";" + specific_length + ";" + str(specific_tm) + ";" + str(specific_gc) + ";" + str(specific_any) + ";" + str(specific_end) + ";" + \
                        str(pair_penalty) + "\n")

      oligo_number += 1



#=====================================================================================================================#
# function primers_triplets_file_writing() : to create a file containing one or several triplets of primers for a SNP #
#=====================================================================================================================#

def primers_triplets_file_writing(fl_fam_primer3_output_parser_name,fl_hex_primer3_output_parser_name) :

  # file containing the primers doublets for the FAM allele
  fl_fam = open(fl_fam_primer3_output_parser_name, 'r')
  # file containing the triplets of primers
  fl_triplets_output = open("primers_triplets.csv", 'w')

  fl_triplets_output.write("id;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")

  line_fam = fl_fam.readline()

  global nb_of_triplets
  nb_of_triplets = 0

  while 1 :
    line_fam = fl_fam.readline()

    if line_fam == "" :
      # EOF ...
      break

    fam_lst_field = line_fam.split('OLIGO_')

    fam_lst_field_02 = fam_lst_field[1].split(";")

    # file containing the primers doublets for the HEX allele
    fl_hex = open(fl_hex_primer3_output_parser_name, 'r')

    line_hex = fl_hex.readline()

    while 1 :
      line_hex = fl_hex.readline()

      if line_hex == "" :
        # EOF ...
        fl_hex.close()
        break

      hex_lst_field = line_hex.split('OLIGO_')

      hex_lst_field_02 = hex_lst_field[1].split(";")

      if fam_lst_field_02[1] == hex_lst_field_02[1] :  # and nb_of_triplets < 5 : # the number of triplets is limited to 5 for each SNP

        score = float(chop(fam_lst_field_02[15])) + float(chop(hex_lst_field_02[15]))

        if score >= 0. and score < 5. :
          category = 1
        elif score >= 5. and score < 10. :
          category = 2
        elif score >= 10. and score < 15. :
          category = 3
        elif score >= 15. and score < 20. :
          category = 4
        elif score >= 20. :
          category = 5

        fl_triplets_output.write(fam_lst_field[0] + fam_lst_field_02[1] + ";" + fam_lst_field_02[2] + ";" + fam_lst_field_02[3] + ";" + fam_lst_field_02[4] + ";" + \
                                                    fam_lst_field_02[5] + ";" + fam_lst_field_02[6] + ";" + fam_lst_field_02[7] + ";" + \
                                                    fam_lst_field_02[8] + ";" + fam_lst_field_02[9] + ";" + fam_lst_field_02[10] + ";" + fam_lst_field_02[11] + ";" + \
                                                    fam_lst_field_02[12] + ";" + fam_lst_field_02[13] + ";" + fam_lst_field_02[14] + ";" + chop(fam_lst_field_02[15]) + ";" + \
                                                    hex_lst_field_02[8] + ";" + hex_lst_field_02[9] + ";" + hex_lst_field_02[10] + ";" + hex_lst_field_02[11] + ";" + \
                                                    hex_lst_field_02[12] + ";" + hex_lst_field_02[13] + ";" + hex_lst_field_02[14] + ";" + chop(hex_lst_field_02[15]) + ";" + \
                                                    str(score) + ";" + str(category) + "\n")

  if fl_fam :
    fl_fam.close()

  if fl_triplets_output :
    fl_triplets_output.close()

  # triplets of primers are sorted by increasing score (the best triplets having the lower score)
  reader = csv.reader(open('primers_triplets.csv'),delimiter=';')

  row = reader.next() # first line of the csv file
  list_reader = []
  for row in reader : # from the second line to the end of the csv file
    list_reader.append(row)

  list_reader.sort(key=lambda v: float(v[24])) # sort by increasing score

  fl_triplets_output_sorted = open("primers_triplets_sorted.csv",'w')
  if len(species_db_name) :
    fl_triplets_output_sorted.write("id;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")
  else :
    fl_triplets_output_sorted.write("id;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")
  i = 0
  while i < min(len(list_reader),5) :  # the number of triplets in primers_triplets_sorted.csv is limited to the top 5 for each SNP
    fl_triplets_output_sorted.write(';'.join(list_reader[i]) + "\n")
    if len(species_db_name) :
      fl_triplet_fasta = open("triplet_" + str(i+1) + ".fasta",'w')
      fl_triplet_fasta.write(">" + list_reader[i][0] + " common primer\n")
      fl_triplet_fasta.write(list_reader[i][1] + "\n")
      fl_triplet_fasta.write(">" + list_reader[i][0] + " FAM specific primer\n")
      fl_triplet_fasta.write(list_reader[i][8] + "\n")
      fl_triplet_fasta.write(">" + list_reader[i][0] + " HEX specific primer\n")
      fl_triplet_fasta.write(list_reader[i][16] + "\n")
      if fl_triplet_fasta :
        fl_triplet_fasta.close()
    i += 1
  nb_of_triplets = i
  if fl_triplets_output_sorted :
    fl_triplets_output_sorted.close()

  # Adding a column "triplet_nb" in file primers_triplets_list_sorted.txt
  os.rename("primers_triplets_sorted.csv", "primers_triplets_sorted_tmp.csv")
  fl_triplets_output_sorted_tmp = open("primers_triplets_sorted_tmp.csv",'r')
  fl_triplets_output_sorted = open("primers_triplets_sorted.csv",'w')
  if len(species_db_name) :
    fl_triplets_output_sorted.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")
  else :
    fl_triplets_output_sorted.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")
  triplet_nb = 1
  line = fl_triplets_output_sorted_tmp.readline()

  while 1 :
    line = fl_triplets_output_sorted_tmp.readline()
    if line == "" :
        # EOF ...
        break
    lst_field = line.split(";")
    line_fl_triplets_output_sorted = fam_lst_field[0] + "triplet_" + str(triplet_nb) + ";"
    i = 0
    while i < len(lst_field) - 2 :
      line_fl_triplets_output_sorted = line_fl_triplets_output_sorted + lst_field[1 + i] + ";"
      i += 1
    # adding the triplet's class (ie the last element of the line)
    line_fl_triplets_output_sorted = line_fl_triplets_output_sorted + lst_field[len(lst_field) - 1]
    fl_triplets_output_sorted.write(line_fl_triplets_output_sorted)

    triplet_nb += 1

  if fl_triplets_output_sorted_tmp :
    fl_triplets_output_sorted_tmp.close()
    os.remove("primers_triplets_sorted_tmp.csv")

  if fl_triplets_output_sorted :
    fl_triplets_output_sorted.close()



#=======================================================================================================#
# function writing_all_primers_triplets_files() : to write files containing all the primers of triplets #
#=======================================================================================================#

def writing_all_primers_triplets_files() :

  # writing a file containing all the triplets of primers of the forward alleles

  os.chdir(current_working_directory + "SNP")
  fl_forward_primers_triplets = open("forward_primers_triplets.csv",'w')
  if len(species_db_name) :
    fl_forward_primers_triplets.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")
  else :
    fl_forward_primers_triplets.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")

  i = 1

  while i <= nb_of_snp :
    # the SNP is processed
    if snp_processing[i-1] == 1 :
      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/FORWARD")
      fl_primers_triplets_list_sorted = open("primers_triplets_sorted.csv").readlines()[1:]
      fl_forward_primers_triplets.writelines(fl_primers_triplets_list_sorted)
    i += 1

  if fl_forward_primers_triplets :
    fl_forward_primers_triplets.close()

  # writing a file containing all the triplets of primers of the reverse alleles

  os.chdir(current_working_directory + "SNP/")
  fl_reverse_primers_triplets = open("reverse_primers_triplets.csv",'w')
  if len(species_db_name) :
    fl_reverse_primers_triplets.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")
  else :
    fl_reverse_primers_triplets.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")

  i = 1

  while i <= nb_of_snp :
    # the SNP is processed
    if snp_processing[i-1] == 1 :
      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/REVERSE")
      fl_primers_triplets_list_sorted = open("primers_triplets_sorted.csv").readlines()[1:]
      fl_reverse_primers_triplets.writelines(fl_primers_triplets_list_sorted)
    i += 1

  if fl_reverse_primers_triplets :
    fl_reverse_primers_triplets.close()

  # writing a file containing all the triplets of primers for the forward and reverse alleles

  os.chdir(current_working_directory + "SNP")
  fl_forward_reverse_primers_triplets = open("all_primers_triplets.csv",'w')
  if len(species_db_name) :
    fl_forward_reverse_primers_triplets.write("id;FORWARD/REVERSE;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")
  else :
    fl_forward_reverse_primers_triplets.write("id;FORWARD/REVERSE;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class \n")

  fl_forward_primers_triplets = open("forward_primers_triplets.csv")
  fl_forward_primers_triplets.readline()
  while 1 :
    line = fl_forward_primers_triplets.readline()
    if line == "" :
      # EOF ...
      break
    lst_field = line.split(";")
    line_fl_primers_triplets = lst_field[0] + ";FORWARD;" + lst_field[1]
    i = 2
    while i < len(lst_field) :
      line_fl_primers_triplets = line_fl_primers_triplets + ";" + lst_field[i]
      i += 1
    fl_forward_reverse_primers_triplets.write(line_fl_primers_triplets)
  if fl_forward_primers_triplets :
    fl_forward_primers_triplets.close()

  fl_reverse_primers_triplets = open("reverse_primers_triplets.csv")
  fl_reverse_primers_triplets.readline()
  while 1 :
     line = fl_reverse_primers_triplets.readline()
     if line == "" :
       # EOF ...
       break
     lst_field = line.split(";")
     line_fl_primers_triplets = lst_field[0] + ";REVERSE;" + lst_field[1]
     i = 2
     while i < len(lst_field) :
       line_fl_primers_triplets = line_fl_primers_triplets + ";" + lst_field[i]
       i += 1
     fl_forward_reverse_primers_triplets.write(line_fl_primers_triplets)
  if fl_reverse_primers_triplets :
    fl_reverse_primers_triplets.close()

  if fl_forward_reverse_primers_triplets :
    fl_forward_reverse_primers_triplets.close()



#====================================================================================================================#
# function testing_primers_triplets_specificities() : to test the specificities of the triplets of primers for a SNP #
#====================================================================================================================#

def testing_primers_triplets_specificities() :

  global list_triplets_specificities

  i = 1

  while i <= nb_of_snp :

    # the SNP is processed
    if snp_processing[i-1] == 1 :

      # testing triplets of primers of SNP i in the forward direction
      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/FORWARD")

      j = 1

      list_triplets_specificities = []  # specificity of each triplet of primers of the SNP i in the forward direction

      while j <= forward_nb_of_triplets[i-1] and forward_nb_of_triplets[i-1] > 0 :

        # testing the specificity of the triplet j for SNP i in the forward direction
        primers_triplet_specificity(i,j,"FORWARD")

        j += 1

      # adding some informations about the specificity of the triplets of primers in file primers_triplets_sorted.csv
      adding_specificities_informations()

      # testing triplets of primers of SNP i in the reverse direction
      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/REVERSE")

      j = 1

      list_triplets_specificities = []  # specificity of each triplet of primers of the SNP i in the reverse direction

      while j <= reverse_nb_of_triplets[i-1] and reverse_nb_of_triplets[i-1] > 0 :

        # testing the specificity of the triplet j for SNP i in the reverse direction
        primers_triplet_specificity(i,j,"REVERSE")

        j += 1

      # adding some informations about the specificity of the triplets of primers in file primers_triplets_sorted.csv
      adding_specificities_informations()

    i += 1

  # writing primers_triplets_specificities.csv file in the SNP directory
  writing_primers_triplets_specificities_file()



#==========================================================================================================#
# function primers_triplet_specificity(i,j,direction) : to test the specificity of the triplet j for SNP i #
#==========================================================================================================#

def primers_triplet_specificity(i,j,direction) :

  os.system("blastall -p blastn -d " + species_db_name + " -i triplet_" + str(j) + ".fasta -F F -e 1.0 -b 250 -v 250 -o triplet_" + str(j) + "_vs_base.bltn")
  os.system("/home/gjb99/GBS/tools/KaSpaR/filtre_blastall_ncbi -nb 100 -all < triplet_" + str(j) + "_vs_base.bltn > triplet_" + str(j) + "_vs_base.bltn.filtre")

  fl_triplet_filtre = open("triplet_" + str(j) + "_vs_base.bltn.filtre")

  etat = 0
  primer_type = ""
  list_common = []
  list_fam = []
  list_hex = []

  while 1 :
    line = fl_triplet_filtre.readline()

    if line == "" :
      # EOF ...
      break

    if etat == 0 :

      if line.find("common") > 0 :

        line = fl_triplet_filtre.readline()
        common_scaffold = line.split(' ')[1]
        primer_type = "common"

      elif line.find("FAM") > 0 :

        line = fl_triplet_filtre.readline()
        fam_scaffold = line.split(' ')[1]
        primer_type = "fam"

      elif line.find("HEX") > 0 :

        line = fl_triplet_filtre.readline()
        hex_scaffold = line.split(' ')[1]
        primer_type = "hex"

      etat = 1

    elif etat == 1 and primer_type == "common" :

      if line == '\n' :

        etat = 0

      else :

        lst_field = line.split(' ')

        common_subject_start = lst_field[4]
        common_subject_end = lst_field[5]
        common_subject_length = lst_field[6]
        common_subject_sens = lst_field[7]
        common_query_start = lst_field[11]
        common_query_end = lst_field[12]
        common_query_length = lst_field[13]
        common_query_sens = lst_field[14]

        if int(common_query_start) == 1 and int(common_query_end) == int(common_query_length) :

          dic_common = {}
          dic_common['scaffold'] = common_scaffold
          dic_common['subject_start'] = common_subject_start
          dic_common['subject_end'] = common_subject_end
          dic_common['subject_length'] = common_subject_length
          dic_common['subject_sens'] = common_subject_sens
          dic_common['query_start'] = common_query_start
          dic_common['query_end'] = common_query_end
          dic_common['query_length'] = common_query_length
          dic_common['query_sens'] = common_query_sens

          list_common.append(dic_common)

    elif etat == 1 and primer_type == "fam" :

      if line == '\n' :

        etat = 0

      else :

        lst_field = line.split(' ')

        fam_subject_start = lst_field[4]
        fam_subject_end = lst_field[5]
        fam_subject_length = lst_field[6]
        fam_subject_sens = lst_field[7]
        fam_query_start = lst_field[11]
        fam_query_end = lst_field[12]
        fam_query_length = lst_field[13]
        fam_query_sens = lst_field[14]

        if int(fam_query_start) == 1 and \
           (int(fam_query_end) == int(fam_query_length) or int(fam_query_end) == int(fam_query_length) -1) :

          dic_fam = {}
          dic_fam['scaffold'] = fam_scaffold
          dic_fam['subject_start'] = fam_subject_start
          dic_fam['subject_end'] = fam_subject_end
          dic_fam['subject_length'] = fam_subject_length
          dic_fam['subject_sens'] = fam_subject_sens
          dic_fam['query_start'] = fam_query_start
          dic_fam['query_end'] = fam_query_end
          dic_fam['query_length'] = fam_query_length
          dic_fam['query_sens'] = fam_query_sens

          list_fam.append(dic_fam)

    elif etat == 1 and primer_type == "hex" :

      if line == '\n' :

        etat = 0

      else :

        lst_field = line.split(' ')

        hex_subject_start = lst_field[4]
        hex_subject_end = lst_field[5]
        hex_subject_length = lst_field[6]
        hex_subject_sens = lst_field[7]
        hex_query_start = lst_field[11]
        hex_query_end = lst_field[12]
        hex_query_length = lst_field[13]
        hex_query_sens = lst_field[14]

        if int(hex_query_start) == 1 and \
           (int(hex_query_end) == int(hex_query_length) or int(hex_query_end) == int(hex_query_length) -1) :

          dic_hex = {}
          dic_hex['scaffold'] = hex_scaffold
          dic_hex['subject_start'] = hex_subject_start
          dic_hex['subject_end'] = hex_subject_end
          dic_hex['subject_length'] = hex_subject_length
          dic_hex['subject_sens'] = hex_subject_sens
          dic_hex['query_start'] = hex_query_start
          dic_hex['query_end'] = hex_query_end
          dic_hex['query_length'] = hex_query_length
          dic_hex['query_sens'] = hex_query_sens

          list_hex.append(dic_hex)

  if fl_triplet_filtre :
    fl_triplet_filtre.close()

  # writing file triplet_j_specificity.csv in SNP/SNP_i/FORWARD and SNP/SNP_i/REVERSE (i = SNP's index, j = triplet's index)
  fl_triplet_specificity = open("triplet_" + str(j) + "_specificity.csv",'w')
  fl_triplet_specificity.write("id;triplet;" + \
                               "scaffold;scaffold's length;SNP's position;" + \
                               "common primer's global start position;common primer's global end position;common primer's global sens;" + \
                               "common primer's local start position;common primer's local end position;common primer's local length;common primer's local sens;" + \
                               "FAM primer's global start position;FAM primer's global end position;FAM primer's global sens;" + \
                               "FAM primer's local start position;FAM primer's local end position;FAM primer's local length;FAM primer's local sens;" + \
                               "HEX primer's global start position;HEX primer's global end position;HEX primer's global sens;" + \
                               "HEX primer's local start position;HEX primer's local end position;HEX primer's local length;HEX primer's local sens \n")

  nb_of_specificities = 0

  if list_common != [] and list_fam != [] and list_hex != [] :
    for dic_c in list_common :
      for dic_f in list_fam :
        for dic_h in list_hex :
          if dic_f['scaffold'] == dic_c['scaffold'] and dic_h['scaffold'] == dic_c['scaffold'] and \
             dic_f['subject_sens'] != dic_c['subject_sens'] and dic_h['subject_sens'] != dic_c['subject_sens'] and \
             abs(int(dic_f['subject_end'])-int(dic_c['subject_end'])) < 1000 and abs(int(dic_h['subject_end'])-int(dic_c['subject_end'])) < 1000 and \
             abs(int(dic_h['subject_end'])-int(dic_f['subject_end'])) == 1 :

            nb_of_specificities += 1

  if nb_of_specificities == 1 :

    for dic_c in list_common :
      for dic_f in list_fam :
        for dic_h in list_hex :
          if dic_f['scaffold'] == dic_c['scaffold'] and dic_h['scaffold'] == dic_c['scaffold'] and \
             dic_f['subject_sens'] != dic_c['subject_sens'] and dic_h['subject_sens'] != dic_c['subject_sens'] and \
             abs(int(dic_f['subject_end'])-int(dic_c['subject_end'])) < 1000 and abs(int(dic_h['subject_end'])-int(dic_c['subject_end'])) < 1000 and \
             abs(int(dic_h['subject_end'])-int(dic_f['subject_end'])) == 1 :

            # scaffold
            scaffold = dic_c['scaffold']

            # SNP's position
            if direction == "FORWARD" :
              # forward direction
              snp_position = str(max(int(dic_f['subject_end']),int(dic_h['subject_end'])))
            else :
              # reverse direction
              snp_position = str(min(int(dic_f['subject_end']),int(dic_h['subject_end'])))

            fl_triplet_specificity.write(lst_snp_id[i-1] + ";triplet_" + str(j) + ";" + \
                                         dic_c['scaffold'] + ";" + \
                                         dic_c['subject_length'] + ";" + \
                                         snp_position + ";" + \
                                         dic_c['subject_start'] + ";" + dic_c['subject_end'] + ";" + dic_c['subject_sens'] + ";" + \
                                         dic_c['query_start'] + ";" + dic_c['query_end'] + ";" + dic_c['query_length'] + ";" + dic_c['query_sens'] + ";" + \
                                         dic_f['subject_start'] + ";" + dic_f['subject_end'] + ";" + dic_f['subject_sens'] + ";" + \
                                         dic_f['query_start'] + ";" + dic_f['query_end'] + ";" + dic_f['query_length'] + ";" + dic_f['query_sens'] + ";" + \
                                         dic_h['subject_start'] + ";" + dic_h['subject_end'] + ";" + dic_h['subject_sens'] + ";" + \
                                         dic_h['query_start'] + ";" + dic_h['query_end'] + ";" + dic_h['query_length'] + ";" + dic_h['query_sens'] + "\n")

  if fl_triplet_specificity :
    fl_triplet_specificity.close()

  dic_triplets_specificities = {}
  if nb_of_specificities == 1 and len(list_common) == 1 and len(list_fam) == 1 and len(list_hex) == 1 :
    dic_triplets_specificities['specific_or_not'] = 'specific'
    dic_triplets_specificities['scaffold'] = scaffold
    dic_triplets_specificities['SNP_position'] = snp_position
  elif nb_of_specificities == 1 and (len(list_common) > 1  or len(list_fam) > 1 or len(list_hex) > 1) :
    dic_triplets_specificities['specific_or_not'] = 'specific ?'
    dic_triplets_specificities['scaffold'] = scaffold
    dic_triplets_specificities['SNP_position'] = snp_position
  else :
    dic_triplets_specificities['specific_or_not'] = 'not specific'
    dic_triplets_specificities['scaffold'] = 'no scaffold'
    dic_triplets_specificities['SNP_position'] = 'no SNP\'s position'

  list_triplets_specificities.append(dic_triplets_specificities)



#==============================================================================================================================================================#
# function adding_specificities_informations() : to add some informations about the specificity of the triplets of primers in file primers_triplets_sorted.csv #
#==============================================================================================================================================================#

def adding_specificities_informations() :

  j = 1

  os.rename("primers_triplets_sorted.csv","primers_triplets_sorted_tmp.csv")
  fl_primers_triplets_sorted_tmp = open("primers_triplets_sorted_tmp.csv")

  fl_primers_triplets_sorted = open("primers_triplets_sorted.csv",'w')

  fl_primers_triplets_sorted.write("id;triplet;common primer;start;len;tm;%gc;any;3prime;FAM specific primer;start;len;tm;%gc;any;3prime;FAM pair's score;HEX specific primer;start;len;tm;%gc;any;3prime;HEX pair's score;triplet's score;triplet's class;specific or not;scaffold;SNP's position \n")

  line = fl_primers_triplets_sorted_tmp.readline()

  while 1 :
    line = fl_primers_triplets_sorted_tmp.readline()
    if line == "" :
      # EOF ...
      break
    fl_primers_triplets_sorted.write(line[:-1] + ";" + list_triplets_specificities[j-1]['specific_or_not'] + ";" + list_triplets_specificities[j-1]['scaffold'] + ";" + \
                                                       list_triplets_specificities[j-1]['SNP_position'] + "\n")
    j += 1

  if fl_primers_triplets_sorted_tmp :
    fl_primers_triplets_sorted_tmp.close()
    os.remove("primers_triplets_sorted_tmp.csv")

  if fl_primers_triplets_sorted :
    fl_primers_triplets_sorted.close()

  # These informations are consequently added in files :
  #  - SNP/forward_primers_triplets.csv
  #  - SNP/reverse_primers_triplets.csv
  #  - SNP/all_primers_triplets.csv
  # (see writing_all_primers_triplets_files())



#========================================================================================================================================================#
# function writing_primers_triplets_specificities_file() : to write file SNP/primers_triplets_specificities.csv including all specificities informations #
#========================================================================================================================================================#

def writing_primers_triplets_specificities_file() :

  os.chdir(current_working_directory + "SNP")

  fl_primers_triplets_specificities = open("primers_triplets_specificities.csv",'w')

  fl_primers_triplets_specificities.write("id;FORWARD/REVERSE;triplet;" + \
                                          "scaffold;scaffold's length;SNP's position;" + \
                                          "common primer's global start position;common primer's global end position;common primer's global sens;" + \
                                          "common primer's local start position;common primer's local end position;common primer's local length;common primer's local sens;" + \
                                          "FAM primer's global start position;FAM primer's global end position;FAM primer's global sens;" + \
                                          "FAM primer's local start position;FAM primer's local end position;FAM primer's local length;FAM primer's local sens;" + \
                                          "HEX primer's global start position;HEX primer's global end position;HEX primer's global sens;" + \
                                          "HEX primer's local start position;HEX primer's local end position;HEX primer's local length;HEX primer's local sens \n")

  i = 1

  # forward direction
  while i <= nb_of_snp :

    # the SNP is processed
    if snp_processing[i-1] == 1 :

      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/FORWARD")

      j = 1

      while j <= forward_nb_of_triplets[i-1] and forward_nb_of_triplets[i-1] > 0 :

        fl_triplet_specificity = open("triplet_" + str(j) + "_specificity.csv",'r')
        fl_triplet_specificity.readline()
        while 1 :
          line = fl_triplet_specificity.readline()
          if line == "" :
            # EOF ...
            break
          lst_field = line.split(";")
          line_fl_primers_triplets_specificities = lst_field[0] + ";FORWARD"
          k = 1
          while k < len(lst_field) :
            line_fl_primers_triplets_specificities = line_fl_primers_triplets_specificities + ";" + lst_field[k]
            k += 1
          fl_primers_triplets_specificities.write(line_fl_primers_triplets_specificities)

        if fl_triplet_specificity :
          fl_triplet_specificity.close()

        j += 1

    i += 1

  i = 1

  # reverse direction
  while i <= nb_of_snp :

    # the SNP is processed
    if snp_processing[i-1] == 1 :

      os.chdir(current_working_directory + "SNP/SNP_" + str(i) + "/REVERSE")

      j = 1

      while j <= reverse_nb_of_triplets[i-1] and reverse_nb_of_triplets[i-1] > 0 :

        fl_triplet_specificity = open("triplet_" + str(j) + "_specificity.csv",'r')
        fl_triplet_specificity.readline()
        while 1 :
          line = fl_triplet_specificity.readline()
          if line == "" :
            # EOF ...
            break
          lst_field = line.split(";")
          line_fl_primers_triplets_specificities = lst_field[0] + ";REVERSE"
          k = 1
          while k < len(lst_field) :
            line_fl_primers_triplets_specificities = line_fl_primers_triplets_specificities + ";" + lst_field[k]
            k += 1
          fl_primers_triplets_specificities.write(line_fl_primers_triplets_specificities)

        if fl_triplet_specificity :
          fl_triplet_specificity.close()

        j += 1

    i += 1

  if fl_primers_triplets_specificities :
    fl_primers_triplets_specificities.close()



#=======================================================================#
# function print_msg(msg) : to display information during the execution #
#=======================================================================#

def print_msg(msg = "\n") :
  if type(msg) == type('a') :
    if msg[len(msg) - 1] != '\n' :
      msg += '\n'
    if fl_err :
      fl_err.write(msg)
      fl_err.flush()
    else :
      sys.stderr.write(msg)
      sys.stderr.flush()

  return



#===========================================================================#
# Main part of this program                                                 #
#===========================================================================#

if __name__ == "__main__" :

  # Defining global variables
  #==========================

  debug                          = 0                   # debug mode [0|1]
  message                        = ""
  nb_arg                         = len(sys.argv)       # count of arguments
  fl_input                       = 0                   # file containing a list of SNP with their backward and forward sequences (-in parameter)
  fl_input_name  = ""
  fl_pcr_parameters_default      = 0                   # file containing the default PCR parameters
  fl_pcr_parameters_default_name = "/usr/local/prive/bi/gafl/pcr_parameters_default.txt"
  fl_pcr_parameters              = 0                   # file containing the PCR parameters
  fl_pcr_parameters_name         = ""
  species_db_name                = ""                  # path to the data base of the species used to check the specificities of the triplets of primers
  fl_err                         = 0                   # error file (-err parameter)
  fl_err_name                    = ""
  nb_of_snp                      = 0                   # number of SNPs
  lst_snp_id                     = []                  # list containing the SNPs identifiers
  current_working_directory      = os.getcwd() + "/"   # current working directory
  flanking_sequence_min_size     = 50                  # minimum length of the flanking sequences
  checked_zone_length            = 10                  # no IUPAC code among the 10 nucleotides before and after the SNP in a sequence
  snp_processing                 = []                  # list indicating if a SNP is processed or not
  forward_nb_of_triplets         = []                  # list indicating the number of triplets of primers for each SNP in the forward direction
  reverse_nb_of_triplets         = []                  # list indicating the number of triplets of primers for each SNP in the reverse direction

  # Processing  arguments
  #======================

  if (nb_arg < 3 or nb_arg > 25) and (nb_arg != 2 or sys.argv[1] not in ("-h", "-H", "-pcr_parameters_default")) :
    error(3)

  list_pcr_parameters = []

  i = 1
  while i < nb_arg :
    if sys.argv[i] == "-h" :
      help()
    elif sys.argv[i] == "-H" :
      help_detailed()
    elif sys.argv[i] == "-pcr_parameters_default" :
      getting_default_pcr_parameters()
    elif sys.argv[i] == "-in" and i < nb_arg :
      fl_input_name = sys.argv[i + 1]
      i += 1
    elif sys.argv[i] == "-flanking_sequence_min_size" and i < nb_arg :
      flanking_sequence_min_size = int(sys.argv[i + 1])
      i += 1
    elif sys.argv[i] == "-checked_zone_length" and i < nb_arg :
      checked_zone_length = int(sys.argv[i + 1])
      i += 1
    elif sys.argv[i] == "-pcr_parameters_file" and i < nb_arg :
      fl_pcr_parameters_name = sys.argv[i + 1]
      i += 1
    elif sys.argv[i] == "-pcr_parameters" and i < nb_arg :
      i += 1
      while sys.argv[i] != "-h" and sys.argv[i] != "-H" and sys.argv[i] != "-pcr_parameters_default" \
                                and sys.argv[i] != "-in" and sys.argv[i] != "-flanking_sequence_min_size" and sys.argv[i] != "-checked_zone_length" \
                                and sys.argv[i] != "-pcr_parameters_file" and sys.argv[i] != "-species_db" and sys.argv[i] != "-err" and sys.argv[i] != "-debug" :
        list_pcr_parameters.append(sys.argv[i])
        if i + 1 == nb_arg or sys.argv[i+1] == "-h" or sys.argv[i+1] == "-H" or sys.argv[i+1] == "-pcr_parameters_default" \
                           or sys.argv[i+1] == "-in" or sys.argv[i+1] == "-flanking_sequence_min_size" or sys.argv[i+1] == "-checked_zone_length" \
                           or sys.argv[i+1] == "-pcr_parameters_file" or sys.argv[i+1] == "-species_db" or sys.argv[i+1] == "-err" or sys.argv[i+1] == "-debug" :
          break
        else :
          i += 1
    elif sys.argv[i] == "-species_db" and i < nb_arg :
      species_db_name = sys.argv[i + 1]
      i += 1
    elif sys.argv[i] == "-err" and i < nb_arg :
      fl_err_name = sys.argv[i + 1]
      i += 1
    elif sys.argv[i] == "-debug" :
      debug = 1
    else :
      error(3)
    i += 1

  # Checking arguments and opening files
  # ====================================

  # fl_input_name : file to parse
  if len(fl_input_name) :
    try :
      fl_input = open(fl_input_name, 'r')
    except :
      error(6, fl_input_name)
  else :
      error(3, "-in")

  # errors file
  if len(fl_err_name) :
    try :
      fl_err = open(fl_err_name, 'w')
    except :
      error(7, fl_err_name)

  # Displaying an information during the execution
  # ==============================================

  # arguments seem right
  if debug :
    message = "Processing file " + fl_input_name
    print_msg(message)

  # Designing one or several triplets of primers for each snp
  # =========================================================
  primers_triplets_design()

  # Closing files
  # =============
  files_close()

  sys.exit(0)
