#define _CRT_SECURE_NO_WARNINGS

#include <vector>
#include <stdint.h>
#include <assert.h>
typedef int64_t isize;

struct Edge {
    int32_t src;
    int32_t dest;
    int32_t weight;
};

static isize kruskal(const Edge* sorted_edges, isize edges_count, isize vertex_count)
{
    struct Set {
        int32_t parent;
        int32_t rank;
    };
    
    std::vector<Set> disjoint_sets_vec(vertex_count, Set{-1, 0});
    Set* disjoint_sets = disjoint_sets_vec.data();

    isize weight_sum = 0;
    for(isize i = 0; i < edges_count; i++) 
    {
        Edge edge = sorted_edges[i];

        int32_t i1 = edge.src;
        int32_t i2 = edge.dest;
        assert(0 <= i1 && i1 < vertex_count);
        assert(0 <= i2 && i2 < vertex_count);

        while(disjoint_sets[i1].parent != -1)
            i1 = disjoint_sets[i1].parent;
            
        while(disjoint_sets[i2].parent != -1)
            i2 = disjoint_sets[i2].parent;
            
        assert(0 <= i1 && i1 < vertex_count);
        assert(0 <= i2 && i2 < vertex_count);
        if (i1 != i2) 
        {
            Set* s1 = &disjoint_sets[i1];
            Set* s2 = &disjoint_sets[i2];
            if (s1->rank < s2->rank) {
                s1->parent = i2;
            }
            else if (s1->rank > s2->rank) {
                s2->parent = i1;
            }
            else {
                s2->parent = i1;
                s1->rank += 1;
            }
            
            weight_sum += edge.weight;
        }
    }

    return weight_sum;
}

//The necessary filler bellow... 
// Can this be done in a simpler way? Yes, of course. 
// However, the file parsing takes significant amount of time - somewhere between 100% and 25% of runtime
// of the algorithm itself. Thus we need to be at least little efficient. We do not attempt to present
// "optimized" file parsing, instead only do "non-pessimized" - not perform any unnecessary operations.
// Turns out if we just let cpu do its thing its plenty fast and we can parse 250MB graph file in 250ms.

#include <iostream>
#include <utility>
#include <algorithm>
#include <chrono>

static bool string_skip_deci(const char* data, isize size, isize* i, isize* val)
{
    (void) size; //eh I am leaving all my inputs null terminated so its fine
    isize before = *i;
    isize parsed = 0;
    for(;; *i += 1)
    {
        char c = data[*i];
        uint8_t val = (uint8_t) c - (uint8_t) '0';
        if(val > 9)
            break;

        parsed = parsed*10 + val;
    }

    *val = parsed;
    return before != *i;
}

static bool string_skip_char(const char* data, isize size, isize* i, char c)
{
    (void) size;
    if(data[*i] == c) {
        *i += 1;
        return true;
    }
    return false;
}

static isize string_find_char(const char* data, isize size, isize from, char c)
{
    for(isize i = from; i < size; i++)
        if(data[i] == c)
            return i;

    return -1;
}

static int64_t clock_ns()
{
    using namespace std::chrono;
    return duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
}

//Some modified functions from C I have laying around...
static bool read_entire_file(const char* path, void** out_contents, isize* file_size_or_null)
{
    assert(path && out_contents);

    bool state = true;
    FILE* file = fopen(path, "rb");

    if(!file)
    {
        printf("read_entire_file: couldn't open file '%s'\n", path);
        state = false;
    }

    isize file_size = 0;
    if(state)
    {
        state = state && fseek(file, 0, SEEK_END) == 0;
        file_size = ftell(file);
        state = state && fseek(file, 0, SEEK_SET) == 0;

        if(!state) 
            printf("read_entire_file: failed to get file size '%s'\n", path);
    }
    
    char* contents = NULL;
    if(state)
    {
        contents = (char*) malloc(file_size + 1);
        state = contents != NULL;
        if(!state) 
            printf("read_entire_file: out of memory '%s'\n", path);
    }
    
    if(state)
    {
        state = fread(contents, 1, file_size, file) == file_size;
        contents[file_size] = '\0'; //null terminate for easy printing.

        file_size = ftell(file);
        state = state && fseek(file, 0, SEEK_SET) == 0;

        if(!state) 
            printf("read_entire_file: failed read file '%s'\n", path);
    }
    
    if(file)
        fclose(file);

    *out_contents = contents;
    if(file_size_or_null)
        *file_size_or_null = file_size;
    return state;
}

typedef struct String_Buffer_16 {
    char data[16];
} String_Buffer_16;

static String_Buffer_16 format_nanoseconds(int64_t ns)
{
    int64_t sec = (int64_t) 1000*1000*1000;
    int64_t milli = (int64_t) 1000*1000;
    int64_t micro = (int64_t) 1000;

    int64_t abs = ns > 0 ? ns : -ns;
    String_Buffer_16 out = {0};
    if(abs >= sec)
        snprintf(out.data, sizeof out.data, "%.2lfs", (double) ns / (double) sec);
    else if(abs >= milli)
        snprintf(out.data, sizeof out.data, "%.2lfms", (double) ns / (double) milli);
    else if(abs >= micro)
        snprintf(out.data, sizeof out.data, "%.2lfus", (double) ns / (double) micro);
    else
        snprintf(out.data, sizeof out.data, "%llins", (long long) ns);

    return out;
}

int main(int argc, char** argv)
{
    typedef long long ll;
    bool state = argc > 1;

    char* data = NULL;
    isize size = 0;

    //Read file
    int64_t before_read_entire = clock_ns();
    if(state)
        state = read_entire_file(argv[1], (void**) &data, &size);
    int64_t after_read_entire = clock_ns();

    std::vector<Edge> edges;
    isize V = 0;
    isize E = 0;
    
    //parse file
    int64_t before_parse_file = clock_ns();
    if(state)
    {
        isize i = 0;
        isize MAX_WEIGHT = 16384;
        isize MAX_VERTICES = 1ll << 20;
        isize MAX_EDGES = 1ll << 39;

        state = state 
            && string_skip_deci(data, size, &i, &V)
            && string_skip_char(data, size, &i, ' ')
            && string_skip_deci(data, size, &i, &E)
            && string_skip_char(data, size, &i, '\n');

        if(!state)
        {
            isize line_start = 0;
            isize line_end = string_find_char(data, size, line_start, '\n');
            printf("parse_graph_file: ERROR invalid header: expected '[number] [number]\\n' got '%.*s'\n", (int) (line_end - line_start), data + line_start);
        }

        if(V > MAX_VERTICES || E > MAX_EDGES)
        {
            printf("parse_graph_file: ERROR invalid header: number of edges or vertices is past the maximum allowed!"
                " vertices/edges:%lli/%lli max:%lli/%lli", (ll) V, (ll) E, (ll) MAX_VERTICES, (ll) MAX_EDGES);
            state = false;
        }

        if(state)
        {
            edges.reserve((size_t) V);
            for(isize line_i = 1; line_i <= E; line_i++)
            {
                isize line_start = i;

                isize src; (void) src;
                isize dest; (void) dest;
                isize weight; (void) weight;
                
                state = state 
                    && string_skip_deci(data, size, &i, &src)
                    && string_skip_char(data, size, &i, ' ')
                    && string_skip_deci(data, size, &i, &dest)
                    && string_skip_char(data, size, &i, ' ')
                    && string_skip_deci(data, size, &i, &weight)
                    && string_skip_char(data, size, &i, '\n');

                if(!state || src >= V || dest >= V || weight > MAX_WEIGHT)
                {
                    isize line_end = string_find_char(data, size, line_start, '\n');
                    if(!state)
                        printf("parse_graph_file: ERROR invalid syntax on line %lli line:'%.*s'\n", (ll) line_i, (int) (line_end - line_start), data + line_start);
                    else
                        printf("parse_graph_file: ERROR invalid values on line %lli line:'%.*s'\n", (ll) line_i, (int) (line_end - line_start), data + line_start);
                    break;
                }
                else
                {
                    Edge edge = {0};
                    edge.src = (int32_t) src;
                    edge.dest = (int32_t) dest;
                    edge.weight = (int32_t) weight;

                    edges.push_back(edge);
                }
            }
        }
    }
    int64_t after_parse_file = clock_ns();

    //free already parsed file cause why not
    free(data);
    data = NULL;
    size = 0;

    //run algorithm
    int64_t before_kruskal = clock_ns();
    if(state)
    {
        std::sort(edges.begin(), edges.end(), [](Edge e1, Edge e2){ return e1.weight < e2.weight; });
        isize result = kruskal(edges.data(), edges.size(), V);
        printf("%lli\n", (ll)  result);
    }
    int64_t after_kruskal = clock_ns();

    printf("timing read_entire: %s\n", format_nanoseconds(after_read_entire - before_read_entire).data);
    printf("timing parse_file: %s\n", format_nanoseconds(after_parse_file - before_parse_file).data);
    printf("timing kruskal: %s\n", format_nanoseconds(after_kruskal - before_kruskal).data);
    return state == false;
}
