import copy
def print_instance(n):
  print("template<typename F, typename T,")
  # need to generate something like:
  # typename N1, typename U1,
  # typename N2, typename U2,
  # typename N3, typename U3,
  # typename N4, typename U4,
  # typename N5, typename U5
  for i in range(n):
      line = f"typename N{i}, typename U{i}"
      if i < n-1:
          line += ","
      print(line)

  print(">")
  print("static auto _propagate_error(F a_f,")
  # need to generate something like:
  # std::array<T,5>& a_deviations,
  # const uncertain<N1,U1>& a_a1,
  # const uncertain<N2,U2>& a_a2,
  # const uncertain<N3,U3>& a_a3,
  # const uncertain<N4,U4>& a_a4,
  # const uncertain<N5,U5>& a_a5
  print(f"static_vector<T,{n}>& a_deviations,")
  for i in range(n):
      line = f"const uncertain<N{i},U{i}>& a_a{i}"
      if i < n-1:
          line += ","
      print(line)
  print(")")
  print("{")
  # need to generate something like:
  # auto nominal    = f(a_a1.nominal() , a_a2.nominal() , a_a3.nominal() , a_a4.nominal() , a_a5.nominal());
  # a_deviations[0] = f(a_a1.upper()   , a_a2.nominal() , a_a3.nominal() , a_a4.nominal() , a_a5.nominal()) - nominal;
  # a_deviations[1] = f(a_a1.nominal() , a_a2.upper()   , a_a3.nominal() , a_a4.nominal() , a_a5.nominal()) - nominal;
  # a_deviations[2] = f(a_a1.nominal() , a_a2.nominal() , a_a3.upper()   , a_a4.nominal() , a_a5.nominal()) - nominal;
  # a_deviations[3] = f(a_a1.nominal() , a_a2.nominal() , a_a3.nominal() , a_a4.upper()   , a_a5.nominal()) - nominal;
  # a_deviations[4] = f(a_a1.nominal() , a_a2.nominal() , a_a3.nominal() , a_a4.nominal() , a_a5.upper()  ) - nominal;
  function_call_args = [f"a_a{i}.nominal()" for i in range(n) ]
  print('auto nominal = a_f(',', '.join(function_call_args)+');')
  for i in range(n):
      function_call_args[i] = function_call_args[i].replace('nominal()','upper()')
      print(f'a_deviations[{i}]= a_f(',', '.join(function_call_args)+') - nominal;')
      function_call_args[i] = function_call_args[i].replace('upper()','nominal()')
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
