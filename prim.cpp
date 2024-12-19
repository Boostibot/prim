#define _CRT_SECURE_NO_WARNINGS

#include <iostream>

#include <vector>
#include <utility>
#include <queue>


#include <assert.h>
#define ASSERT(x) assert(x)
#define TEST(x, ...) (!(x) ? printf("TEST("#x") failed with: " __VA_ARGS__), abort() : (void) 0)

typedef int64_t isize;

// Helper structure to represent an edge with weight
struct Edge {
    int src;
    int dest;
    int weight;

    bool operator>(const Edge& other) const {
        return weight > other.weight;
    }
};

struct Half_Edge {
    int dest;
    int weight;
};

using WeightedGraph = std::vector<std::vector<Half_Edge>>;

std::vector<Edge> prim(const WeightedGraph& graph, int startVertex) {
    size_t numVertices = graph.size();
    std::vector<bool> visited(numVertices, false);
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    std::vector<Edge> mst;

    visited[startVertex] = true;

    for (const auto& edge : graph[startVertex]) {
        Edge added = {startVertex, edge.dest, edge.weight};
        pq.push(added);
    }
    
    isize max_size = pq.size();
    isize iters = 0;
    while (pq.empty() == false) {
        if(max_size < pq.size())
            max_size = pq.size();
        Edge curr = pq.top();
        pq.pop();
        
        iters += 1;
        if (visited[curr.dest]) 
            continue;

        visited[curr.dest] = true;
        mst.push_back(curr);

        std::vector<Half_Edge> const& neighbours = graph[curr.dest];

        for (const auto& edge : neighbours) {
            if (!visited[edge.dest]) {
                Edge added = {curr.dest, edge.dest, edge.weight};
                pq.push(added);
            }
        }
    }
    
    printf("iters:%lli cap:%lli edges:%lli\n", iters, max_size, (isize) 0);
    return mst;
}

#define WG_ADD_EDGE_SKIP_RESIZE_CHECK 1
#define WG_ADD_EDGE_SKIP_DUPLICIT_CHECK 2

bool wg_add_edge(WeightedGraph* graph, Edge e, int flags = 0)
{
    if((flags & WG_ADD_EDGE_SKIP_RESIZE_CHECK) == 0) {
        int max = std::max(e.src, e.dest);
        if(graph->size() <= max)
            graph->resize(max + 1);
    }
    
    if((flags & WG_ADD_EDGE_SKIP_DUPLICIT_CHECK) == 0) {
        std::vector<Half_Edge>* src_adjacencies = &(*graph)[e.src];
        for(size_t i = 0; i < src_adjacencies->size(); i++)
            if((*src_adjacencies)[i].dest == e.dest)
                return false;
    }

    {
        std::vector<Half_Edge>* src_adjacencies = &(*graph)[e.src];
        Half_Edge entry = {e.dest, e.weight};
        src_adjacencies->push_back(entry);
    }

    {
        std::vector<Half_Edge>* dest_adjacencies = &(*graph)[e.dest];
        Half_Edge entry = {e.src, e.weight};
        dest_adjacencies->push_back(entry);
    }
    return true;
}



typedef enum {
    READ_ENTIRE_FILE_OK = 0,
    READ_ENTIRE_ERROR_OPEN,
    READ_ENTIRE_ERROR_SEEK,
    READ_ENTIRE_ERROR_OUT_OF_MEM,
    READ_ENTIRE_ERROR_READ_ERROR,
} Read_Entire_File_Error;

bool read_entire_file(const char* path, void** out_contents, isize* file_size_or_null)
{
    ASSERT(path && out_contents);

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

isize string_parse_dec_uint_naive(const char* data, isize size, isize from, isize* end)
{
    (void) size;
    isize parsed = 0;
    isize i = from;
    for(;; i++)
    {
        char c = data[i];
        uint8_t val = (uint8_t) c - (uint8_t) '0';
        if(val > 9)
            break;

        parsed = parsed*10 + val;
    }

    *end = i;
    return parsed;
}

__forceinline
isize string_parse_dec_uint(const char* data, isize size, isize from, isize* end)
{
    (void) size;
    isize parsed = 0;
    isize till = from;
    for(;; till++)
    {
        char c = data[till];
        uint8_t val = (uint8_t) c - (uint8_t) '0';
        if(val > 9)
            break;
    }

    if(till == from + 6)
    {
        parsed = (isize) 0
            + (data[from + 0] - '0') * 100000
            + (data[from + 1] - '0') * 10000
            + (data[from + 2] - '0') * 1000
            + (data[from + 3] - '0') * 100
            + (data[from + 4] - '0') * 10
            + (data[from + 5] - '0');
        
        *end = from + 6;
        return parsed;
    }
    else
    {
        isize i = from;
        for(; i <= till-4; i += 4)
        {
            parsed = parsed * 10000
                + (data[i + 0] - '0') * 1000
                + (data[i + 1] - '0') * 100
                + (data[i + 2] - '0') * 10
                + (data[i + 3] - '0');
        }

        for(; i < till; i++)
        {
            uint8_t val = (uint8_t) data[i] - (uint8_t) '0';
            parsed = parsed*10 + val;
        }

        *end = i;
        return parsed;
    }
}

__forceinline
bool string_skip_dec_uint(const char* data, isize size, isize* i, isize* val)
{
    isize after = 0;
    *val = string_parse_dec_uint_naive(data, size, *i, &after);
    if(after == *i)
        return false;

    *i = after;
    return true;
}

bool string_skip_char(const char* data, isize size, isize* i, char c)
{
    //if(*i >= size)
        //return false;

    if(data[*i] == c) {
        *i += 1;
        return true;
    }
    return false;
}

isize string_find_char(const char* data, isize size, isize from, char c)
{
    for(isize i = from; i < size; i++)
        if(data[i] == c)
            return i;

    return -1;
}

struct Graph_Vertex {
    uint32_t edges_i;
    uint32_t edges_count; 
    uint32_t edges_capacity; 
};

struct Weighted_Graph2 {
    Graph_Vertex* vertices;
    isize vertices_count;

    Half_Edge* half_edges;
    isize half_edges_count;
};

bool parse_graph_file_header(const char* data, isize size, isize* out_i, isize* V, isize* E)
{
    bool state = true;
    
    isize line_start = *out_i;
    state = state && string_skip_dec_uint(data, size, out_i, V);
    state = state && string_skip_char(data, size, out_i, ' ');
    state = state && string_skip_dec_uint(data, size, out_i, E);
    state = state && string_skip_char(data, size, out_i, '\n');

    //todo check V, E bounds
    if(!state)
    {
        isize line_end = string_find_char(data, size, line_start, '\n');
        printf("parse_graph_file: ERROR invalid header: expected '[number] [number]\\n' got '%.*s'\n", (int) (line_end - line_start), data + line_start);
    }

    return state;
}

bool parse_graph_file_line(const char* data, isize size, isize* out_i, isize* src, isize* dest, isize* weight, isize V, isize E, isize line_i)
{
    typedef long long ll;
    isize MAX_WEIGHT = 16384;
    isize line_start = *out_i;

    uint8_t state1 = string_skip_dec_uint(data, size, out_i, src);
    uint8_t state2 = string_skip_char(data, size, out_i, ' ');
    uint8_t state3 = string_skip_dec_uint(data, size, out_i, dest);
    uint8_t state4 = string_skip_char(data, size, out_i, ' ');
    uint8_t state5 = string_skip_dec_uint(data, size, out_i, weight);
    uint8_t state6 = string_skip_char(data, size, out_i, '\n');

    bool out = true;
    bool parse_error = (state1 & state2 & state3 & state4 & state5 & state6) == 0;
    if(parse_error || *src >= V || *dest >= V || *weight > MAX_WEIGHT)
    {
        isize line_end = string_find_char(data, size, line_start, '\n');
        if(parse_error)
            printf("parse_graph_file: ERROR invalid syntax on line %lli line:'%.*s'\n", (ll) line_i, (int) (line_end - line_start), data + line_start);
        else
            printf("parse_graph_file: ERROR invalid values on line %lli line:'%.*s'\n", (ll) line_i, (int) (line_end - line_start), data + line_start);

        out = false;
    }

    return out;
}

bool parse_graph_file(WeightedGraph* graph, const char* data, isize size)
{
    isize i = 0;
    bool state = true;
    
    isize V = 0;
    isize E = 0;
    state = state && parse_graph_file_header(data, size, &i, &V, &E);

    //todo check V, E bounds
    if(!state)
    {
        isize line_start = 0;
        isize line_end = string_find_char(data, size, line_start, '\n');
        printf("parse_graph_file: ERROR invalid header: expected '[number] [number]\\n' got '%.*s'\n", (int) (line_end - line_start), data + line_start);
    }

    if(state)
    {
        graph->resize((size_t) V);
        isize edges_count = 0;
        isize MAX_WEIGHT = 16384;
        for(isize line_i = 1; line_i <= E; line_i++)
        {
            isize line_start = i;

            isize src; (void) src;
            isize dest; (void) dest;
            isize weight; (void) weight;

            if(parse_graph_file_line(data, size, &i, &src, &dest, &weight, V, E, line_i) == false)
            {
                state = false;
                break;
            }

            Edge edge = {0};
            edge.src = (int) src;
            edge.dest = (int) dest;
            edge.weight = (int) weight;

            wg_add_edge(graph, edge, WG_ADD_EDGE_SKIP_RESIZE_CHECK);
        }
    }

    return true;
}

bool parse_graph_file2(Weighted_Graph2* wg, const char* data, isize size)
{
    isize i = 0;
    bool state = true;
    
    isize V = 0;
    isize E = 0;
    state = state && parse_graph_file_header(data, size, &i, &V, &E);

    //todo check V, E bounds
    if(!state)
    {
        isize line_start = 0;
        isize line_end = string_find_char(data, size, line_start, '\n');
        printf("parse_graph_file: ERROR invalid header: expected '[number] [number]\\n' got '%.*s'\n", (int) (line_end - line_start), data + line_start);
    }

    Edge* edges = NULL;
    Half_Edge* half_edges = NULL;
    Graph_Vertex* vertices = NULL;
    if(state)
    {
        edges       = (Edge*)           malloc(sizeof(Edge)*E);
        half_edges  = (Half_Edge*)      malloc(sizeof(Half_Edge)*E*2);
        vertices    = (Graph_Vertex*)   calloc(sizeof(Graph_Vertex),V);
        state = edges && half_edges && vertices;
        if(state == false)
        {
            isize total_B = 0
                + sizeof(Edge)*E
                + sizeof(Half_Edge)*E*2
                + sizeof(Graph_Vertex)+V;
        
            printf("parse_graph_file: ERROR allocation of %.2lfK vertices %.2lfK edges failed requires total %.2lfMB\n", 
                (double)V/1000, (double)E/1000, (double) total_B/1024/1024);
                
            free(edges);
            free(half_edges);
            free(vertices);
        }
    }

    if(state)
    {
        isize edges_count = 0;
        isize MAX_WEIGHT = 16384;
        for(isize line_i = 1; line_i <= E; line_i++)
        {
            isize line_start = i;

            isize src; (void) src;
            isize dest; (void) dest;
            isize weight; (void) weight;

            if(parse_graph_file_line(data, size, &i, &src, &dest, &weight, V, E, line_i) == false)
            {
                state = false;
                break;
            }

            Edge edge = {0};
            edge.src = (int) src;
            edge.dest = (int) dest;
            edge.weight = (int) weight;

            ASSERT(edges_count < E);
            edges[edges_count ++] = edge;

            ASSERT(src < V);
            ASSERT(dest < V);
            vertices[src].edges_capacity += 1;
            vertices[dest].edges_capacity += 1;
        }
    }

    if(state)
    {
        //assign mem
        uint32_t from = 0;
        for(isize k = 0; k < V; k++)
        {
            vertices[k].edges_i = from;
            from += vertices[k].edges_capacity;
        }

        //fill vertices with half edges
        for(isize k = 0; k < E; k++)
        {
            Edge edge = edges[k];

            Graph_Vertex* dest_v = &vertices[edge.dest];
            Graph_Vertex* src_v = &vertices[edge.src];

            Half_Edge* dest_he = half_edges + dest_v->edges_i + dest_v->edges_count;
            Half_Edge* src_he = half_edges + src_v->edges_i + src_v->edges_count;
            
            dest_v->edges_count += 1;
            src_v->edges_count += 1;

            dest_he->dest = edge.src;
            dest_he->weight = edge.weight;

            src_he->dest = edge.dest;
            src_he->weight = edge.weight;

            ASSERT(dest_v->edges_count <= dest_v->edges_capacity);
            ASSERT(src_v->edges_count <= src_v->edges_capacity);
        }

        wg->half_edges = half_edges;
        wg->half_edges_count = 2*E;
        wg->vertices = vertices;
        wg->vertices_count = V;
        free(edges);
    }

    return true;
}


__forceinline
void prim2(Weighted_Graph2 graph, isize startVertex, void (*push_vertex)(uint32_t src, uint32_t dest, uint32_t weight, void* context), void* context) {
    std::vector<bool> visited(graph.vertices_count, false);
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    
    ASSERT(0 <= startVertex && startVertex < graph.vertices_count);
    visited[startVertex] = true;

    //add first edge neighbours
    {
        Graph_Vertex first_vertex = graph.vertices[startVertex];
        Half_Edge* edges = graph.half_edges + first_vertex.edges_i;
        for(uint32_t i = 0; i < first_vertex.edges_count; i++)
        {
            Half_Edge edge = edges[i];
            Edge added = {startVertex, edge.dest, edge.weight};
            pq.push(added);
        }
    }

    isize max_size = pq.size();
    isize iters = 0;
    while (pq.empty() == false) {
        if(max_size < pq.size())
            max_size = pq.size();

        iters += 1;
        Edge curr = pq.top();
        pq.pop();

        if (visited[curr.dest]) 
            continue;

        visited[curr.dest] = true;

        Graph_Vertex curr_vertex = graph.vertices[curr.dest];
        //push_vertex(curr.src, curr.dest, curr.weight, context);

        Half_Edge* neighbours = graph.half_edges + curr_vertex.edges_i;
        for(uint32_t i = 0; i < curr_vertex.edges_count; i++)
        {
            Half_Edge edge = neighbours[i];
            if (!visited[edge.dest]) {
                Edge added = {curr.dest, edge.dest, edge.weight};
                pq.push(added);
            }
        }
    }

    printf("iters:%lli cap:%lli edges:%lli\n", iters, max_size, graph.half_edges_count);
}


struct Weighted_Graph3 {
    std::vector<Edge> edges;
    isize E = 0;
    isize V = 0;
};

bool parse_graph_file3(Weighted_Graph3* graph, const char* data, isize size)
{
    isize i = 0;
    bool state = true;
    
    isize V = 0;
    isize E = 0;
    state = state && parse_graph_file_header(data, size, &i, &V, &E);

    //todo check V, E bounds
    if(!state)
    {
        isize line_start = 0;
        isize line_end = string_find_char(data, size, line_start, '\n');
        printf("parse_graph_file: ERROR invalid header: expected '[number] [number]\\n' got '%.*s'\n", (int) (line_end - line_start), data + line_start);
    }

    if(state)
    {
        graph->E = E;
        graph->V = V;
        graph->edges.reserve((size_t) E);
        isize edges_count = 0;
        isize MAX_WEIGHT = 16384;
        for(isize line_i = 1; line_i <= E; line_i++)
        {
            isize line_start = i;

            isize src; (void) src;
            isize dest; (void) dest;
            isize weight; (void) weight;

            if(parse_graph_file_line(data, size, &i, &src, &dest, &weight, V, E, line_i) == false)
            {
                state = false;
                break;
            }

            Edge edge = {0};
            edge.src = (int) src;
            edge.dest = (int) dest;
            edge.weight = (int) weight;

            graph->edges.push_back(edge);
        }
    }

    return true;
}

isize kruskal(const Edge* sorted_edges, isize edges_count, isize vertex_count)
{
    struct Set {
        int32_t parent;
        int32_t rank;
    };
    
    std::vector<Set> disjoint_sets(vertex_count, Set{-1, 1});

    isize til = vertex_count - 1;
    if(til > edges_count)
        til = edges_count;

    isize weight_sum = 0;
    for(isize i = 0; i < til; i++) 
    {
        Edge edge = sorted_edges[i];

        int32_t i1 = edge.src;
        int32_t i2 = edge.dest;

        while(disjoint_sets[i1].parent != -1)
            i1 = disjoint_sets[i1].parent;
            
        while(disjoint_sets[i2].parent != -1)
            i2 = disjoint_sets[i2].parent;

        if (i1 != i2) 
        {
            ASSERT(0 <= i1 && i1 < vertex_count);
            ASSERT(0 <= i2 && i2 < vertex_count);
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

isize clock_ns()
{
    return clock() * 1'000'000'000/CLOCKS_PER_SEC;
}

double clock_s()
{
    static bool init = false;
    static clock_t start = 0;
    if(init == false)
    {
        start = clock();
        init = true;
    }

    return (double) (clock() - start)/CLOCKS_PER_SEC;
}

int main()
{
    if(1)
    {
        WeightedGraph wg;
        isize entire_file_size = 0;
        char* entire_file = NULL;

        double read_before = clock_s();
        TEST(read_entire_file("testing_big/LiJOeSx5HUE3XJzyygbFsg.txt", (void**) &entire_file, &entire_file_size));
        double read_after = clock_s();

        Weighted_Graph2 wg2 = {0};
        Weighted_Graph3 wg3;

        double parse_before = clock_s();
        //TEST(parse_graph_file(&wg, entire_file, entire_file_size));
        TEST(parse_graph_file3(&wg3, entire_file, entire_file_size));
        double parse_after = clock_s();
        
        double prim_before = clock_s();
        //std::vector<Edge> res = prim(wg, 0);
        double prim_after = clock_s();
        
        double kruskal_before = clock_s();
        std::sort(wg3.edges.begin(), wg3.edges.end(), [](Edge e1, Edge e2){ return e1.weight < e2.weight; });
        kruskal(wg3.edges.data(), wg3.edges.size(), wg3.V);
        //prim2(wg2, 0, NULL, NULL);
        double kruskal_after = clock_s();

        printf("reading: %e ms\n", (read_after - read_before)*1000);
        printf("parsing: %e ms\n", (parse_after - parse_before)*1000);
        printf("prim:     %e ms\n", (prim_after - prim_before)*1000);
        printf("kruskal:  %e ms\n", (kruskal_after - kruskal_before)*1000);

        //std::vector<isize> sizes;
        //isize size_sum = 0;
        //for(size_t i = 0; i < wg.size(); i++)
        //{
        //    sizes.push_back(wg[i].size());
        //    size_sum += wg[i].size();
        //}
        //std::sort(sizes.begin(), sizes.end());
        //printf("avg size: %.2lf\n", (double) (size_sum)/wg.size());
        //printf("median size: %.2lf\n", (double) sizes[sizes.size()/2]);
    }

    if(0)
    {
        printf("started\n");
        WeightedGraph wg;
    
        TEST(wg_add_edge(&wg, Edge{0, 1, 3}));
        TEST(wg_add_edge(&wg, Edge{0, 2, 2}));
        TEST(wg_add_edge(&wg, Edge{0, 3, 1}));
        TEST(wg_add_edge(&wg, Edge{1, 2, 5}));
        TEST(wg_add_edge(&wg, Edge{2, 3, 1}));
    
        printf("done\n");
    
        std::vector<Edge> res = prim(wg, 0);
        (void) res;

        printf("exiting\n");
    }
}
