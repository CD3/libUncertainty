import copy
def print_instance(n):
  print("template<typename F, typename T,")
  # need to generate something like:
  # typename A0,
  # typename A1,
  # typename A2,
  # typename A3,
  # typename A4
  for i in range(n):
      line = f"typename A{i}"
      if i < n-1:
          line += ","
      print(line)

  print(">")
  print("static auto _propagate_error(F a_f,")
  # need to generate something like:
  # std::array<T,5>& a_deviations,
  # const A0& a_a0,
  # const A1& a_a1,
  # const A2& a_a2,
  # const A3& a_a3,
  # const A4& a_a4
  print(f"static_vector<T,{n}>& a_deviations,")
  for i in range(n):
      line = f"const A{i}& a_a{i}"
      if i < n-1:
          line += ","
      print(line)
  print(")")
  print("{")
  # need to generate something like:
  function_call_args = [f"get_nominal(a_a{i})" for i in range(n) ]
  print('auto nominal = a_f(',', '.join(function_call_args)+');')
  for i in range(n):
      function_call_args[i] = function_call_args[i].replace('nominal','upper')
      print(f'if( is_uncertain(a_a{i}) ) '+'{')
      print(f'a_deviations[{i}]= a_f(',', '.join(function_call_args)+') - nominal;')
      print("}else{")
      print(f"a_deviations[{i}] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units")
      print("}")
      function_call_args[i] = function_call_args[i].replace('upper','nominal')
  print("return nominal;")
  print("}")



print("// BEGIN GENERATED CODE")
print("// this code was generated using the generate_basic_error_propagator_propagate_error_templates.py script")
print("// to delete this code, you can run (in vim) :g/^\s*\/\/ BEGIN GENERATED CODE/,/^\s*\/\/ END GENERATED CODE/ d")
print("// to reinsert it, you can run (in vim) :.! python scripts/generate_basic_error_propagator_propagate_error_templates.py")
print("// I'm sorry this is so old-school, but your debugger will thank me...")
for i in range(20):
    print_instance(i+1)
    print()
print("// END GENERATED CODE")
