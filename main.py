import sys

sys.setrecursionlimit(1000000)


class node:
    def __init__(self, kmer, index):
        self.kmer = kmer
        self.in_num = 0
        self.out_num = 0
        self.index = index
        self.cnt = 1
        self.exist = True
        self.succ = list()


hash_str = {}
k = 30
nodes = list()
gene = ["A", "T", "G", "C"]
kmer_cnt = 0


def inverse(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'G':
        return 'C'
    if c == 'C':
        return 'G'
    return 'N'


def dfs(cur, ret, path_used):
    if not cur.exist:
        return ret
    if len(ret) >= 10000:
        return ret
    if not cur.succ:
        return ret
    local_ret = ret
    local_path = path_used
    cur_max = ""
    for next_nod in cur.succ:
        if not local_path.get((cur.index, next_nod)):
            local_ret += nodes[next_nod].kmer[k - 1]
            local_path[(cur.index, next_nod)] = 1
            t =  dfs(nodes[next_nod], local_ret, local_path)
            if len(t) > len(cur_max):
                cur_max = t
    return cur_max


def search_path(start, out_file):  # 搜索所有可能路径，舍弃出现次数低的节点，适用于前三个数据
    ret = ""
    for i in range(k - 1):
        ret += start.kmer[i]
    if not start.kmer.isupper():
        return 0
    ans = dfs(start, ret, {})
    if len(ans) > 8000:
        out_file.writelines([">ok\n",  ans + "\n"]   )

def search_path_greedy(start, out_file):  # 贪心搜索路径，适用于深搜复杂度过高的data4，效果很差，弃用
    cur = start.index
    ret = nodes[cur].kmer[0:k-1]
    path_used = {}
    while len(ret) < 100000 and nodes[cur].out_num > 0:
        ret = ret + nodes[cur].kmer[k - 1]
        cnt = -1
        next_node = 0
        for it in nodes[cur].succ:
            if (not path_used.get((cur, it))) and nodes[it].cnt > cnt:
                cnt = nodes[it].out_num
                next_node = it
        if not next_node:
            break
        path_used[(cur, next_node)] = 1
        cur = next_node
    if len(ret)> 18000:
        out_file.writelines([">ok\n",  ret + "\n"])


def delete_node(cur, min_cnt):
    nodes[cur].exist = False
    for i in nodes[cur].succ:
        nodes[i].in_num -= 1
        if nodes[i].cnt < min_cnt:
            nodes[cur].out_num -= 1
            delete_node(i, min_cnt)


def brush_nodes():  # 处理bubble\tip\erroneous connection
    for i in range(len(nodes)):
        if nodes[i].out_num > 1:
            max_cnt = 0
            for j in nodes[i].succ:
                if nodes[j].cnt > max_cnt:
                    max_cnt = nodes[j].cnt
            for j in nodes[i].succ:
                if nodes[j].cnt <= max_cnt // 10:
                    nodes[i].out_num -= 1
                    delete_node(j, max_cnt // 10)


def main():
    f = open("short_1.fasta", "r")
    lines = f.readlines()
    data = list()
    for line in lines:
        if line[0] == '>':
            continue
        t = line.strip("\n\r")
        length = len(t)
        t1 = ""
        t2 = ""
        t3 = ""
        for i in range(length):
            t1 += t[length - i - 1]
            t2 += inverse(t[i])
            t3 += inverse(t1[i])
        data.append(t)
        data.append(t1)
        data.append(t2)
        data.append(t3)
    f.close()
    f = open("short_2.fasta", "r")
    lines = f.readlines()
    for line in lines:
        if line[0] == '>':
            continue
        t = line.strip("\n\r")
        length = len(t)
        t1 = ""
        t2 = ""
        t3 = ""
        for i in range(length):
            t1 += t[length - i - 1]
            t2 += inverse(t[i])
            t3 += inverse(t1[i])
        data.append(t)
        data.append(t1)
        data.append(t2)
        data.append(t3)
    f.close()
    global kmer_cnt
    kmer_cnt = 0
    global nodes
    for dat in data:
        L = len(dat)
        for j in range(L - k + 1):
            s = dat[j: k + j]
            if not hash_str.get(s):
                hash_str[s] = kmer_cnt
                t = node(kmer=s, index=kmer_cnt)
                kmer_cnt += 1
                nodes.append(t)
            else:
                nodes[hash_str[s]].cnt += 1

    print(len(nodes))
    for nod in nodes:
        tmp = nod.kmer[1: k]
        for i in gene:
            pp = tmp + i
            if hash_str.get(pp):
                tt = hash_str[pp]
                nod.out_num += 1
                nod.succ.append(tt)
                nodes[tt].in_num += 1
    print("finished making graph")
    brush_nodes()
    f = open("result.txt", "w")
    tot = 0
    for nod in nodes:
        if nod.in_num == 0:
            tot += 1
            print(tot)
            search_path(nod, f)

            if tot > 1000: #限制输出规模，减少重复率和文件大小
                break

    return 0


if __name__ == "__main__":
    main()

