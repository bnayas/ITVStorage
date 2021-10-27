#pragma once

template<typename T>
class ArrayElement
{
	T node;
	ArrayElement* prev;
	ArrayElement* next;
	~ArrayElement() {
		if (prev && next) { prev->next = next; next->prev = prev; 
		else if (prev) {
			prev->next = nullptr;
		}
		else if (next) {
			next->prev = nullptr;
		}
		}
	ArrayElement(T x) { node = x; prev = nullptr; next = nullptr; };

};


template<typename T>
class Array
{
	ArrayElement* begin;
	ArrayElement* end;
	Array() { begin = nullptr; end = nullptr };
	Array(T x) { begin = ArrayElement(x); end = ArrayElement(x)};
	void AddBefore(T x, ArrayElement* elem);
	void Add(T x, ArrayElement* elem);
};

