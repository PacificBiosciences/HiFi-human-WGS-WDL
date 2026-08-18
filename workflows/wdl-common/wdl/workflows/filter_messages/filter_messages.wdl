version 1.0

workflow filter_messages {
  meta {
    description: "Flatten message arrays and drop empty-string placeholders."
    category: "Utility"
    outputs: {
      messages: {
        description: "Flattened, non-empty messages"
      }
    }
  }

  parameter_meta {
    message_arrays: {
      description: "Message arrays to flatten and filter"
    }
  }

  input {
    #@ except: DeclarationName
    Array[Array[String]] message_arrays
  }

  scatter (m in flatten(message_arrays)) {
    if (m != "") {
      String non_empty = m
    }
  }

  output {
    Array[String] messages = select_all(non_empty)
  }
}

