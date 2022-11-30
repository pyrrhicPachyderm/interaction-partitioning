#ifndef UTILS_ARRAY_HPP
#define UTILS_ARRAY_HPP

#include <array>
#include <utility>
#include <type_traits>

//Applies func to each element of arr, returning an array of equal length to the input array.
template<typename Func, typename InputType, size_t n> std::array<std::invoke_result_t<Func, InputType>, n> array_map(Func func, std::array<InputType, n> arr) {
	return [] <std::size_t... I> (Func func, std::array<InputType, n> arr, std::index_sequence<I...>) -> std::array<std::invoke_result_t<Func, InputType>, n> {
		return {func(arr[I])...};
	} (func, arr, std::make_index_sequence<n>{});
}

//As make_index_sequence, but returning std::array instead of std::index_sequence.
template<size_t n> constexpr std::array<size_t, n> make_index_array() {
	return [] <std::size_t... I> (std::index_sequence<I...>) -> std::array<size_t, n> {
		return {I...};
	} (std::make_index_sequence<n>{});
}

#endif
