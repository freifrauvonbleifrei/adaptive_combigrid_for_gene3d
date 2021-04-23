from Qes_data import Qes_data

if __name__ == "__main__":
    lvs = [[5,5,9,5,3],[5,5,6,9,3],[5,5,6,5,6]]
    for lv in lvs:
        qes_data = Qes_data(lv)
        qes = qes_data.get_result()
        frame = qes_data.get_frame()
        print(frame)
        print(qes)
        qes_data.set_csv()